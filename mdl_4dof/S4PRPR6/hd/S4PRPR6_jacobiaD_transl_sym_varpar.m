% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRPR6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(6)) * t14, 0, 0; 0, sin(pkin(6)) * t14, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:25
	% EndTime: 2019-12-29 12:09:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->7), mult. (56->16), div. (0->0), fcn. (41->6), ass. (0->9)
	t123 = r_i_i_C(3) + qJ(3);
	t119 = cos(qJ(2));
	t122 = qJD(2) * t119;
	t121 = -r_i_i_C(1) * cos(pkin(7)) + r_i_i_C(2) * sin(pkin(7)) - pkin(2);
	t118 = sin(qJ(2));
	t120 = qJD(3) * t119 + (-t123 * t118 + t121 * t119) * qJD(2);
	t117 = cos(pkin(6));
	t115 = sin(pkin(6));
	t1 = [0, t120 * t117, t117 * t122, 0; 0, t120 * t115, t115 * t122, 0; 0, t118 * qJD(3) + (t121 * t118 + t123 * t119) * qJD(2), qJD(2) * t118, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:30
	% EndTime: 2019-12-29 12:09:30
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (68->23), mult. (129->49), div. (0->0), fcn. (100->7), ass. (0->18)
	t157 = r_i_i_C(3) + pkin(5) + qJ(3);
	t143 = sin(pkin(6));
	t147 = cos(qJ(2));
	t156 = t143 * t147;
	t144 = cos(pkin(6));
	t155 = t144 * t147;
	t154 = qJD(2) * t147;
	t146 = sin(qJ(2));
	t153 = qJD(4) * t146;
	t142 = pkin(7) + qJ(4);
	t140 = sin(t142);
	t141 = cos(t142);
	t152 = r_i_i_C(1) * t140 + r_i_i_C(2) * t141;
	t151 = -r_i_i_C(1) * t141 + r_i_i_C(2) * t140 - cos(pkin(7)) * pkin(3) - pkin(2);
	t150 = t152 * t146;
	t149 = qJD(2) * t150;
	t148 = qJD(3) * t147 + qJD(4) * t150 + (-t157 * t146 + t151 * t147) * qJD(2);
	t1 = [0, t148 * t144, t144 * t154, t144 * t149 + ((-t140 * t143 - t141 * t155) * r_i_i_C(1) + (t140 * t155 - t141 * t143) * r_i_i_C(2)) * qJD(4); 0, t148 * t143, t143 * t154, t143 * t149 + ((t140 * t144 - t141 * t156) * r_i_i_C(1) + (t140 * t156 + t141 * t144) * r_i_i_C(2)) * qJD(4); 0, t146 * qJD(3) - t152 * t147 * qJD(4) + (t151 * t146 + t157 * t147) * qJD(2), qJD(2) * t146, (t140 * t153 - t141 * t154) * r_i_i_C(2) + (-t140 * t154 - t141 * t153) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end