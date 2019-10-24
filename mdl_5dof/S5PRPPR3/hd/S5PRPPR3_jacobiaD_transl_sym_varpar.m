% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPPR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->5), mult. (24->10), div. (0->0), fcn. (15->6), ass. (0->5)
	t15 = qJ(2) + pkin(8);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = qJD(2) * (-cos(qJ(2)) * pkin(2) - r_i_i_C(1) * t14 + r_i_i_C(2) * t13);
	t1 = [0, cos(pkin(7)) * t19, 0, 0, 0; 0, sin(pkin(7)) * t19, 0, 0, 0; 0, (-sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:53
	% EndTime: 2019-10-24 10:22:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->9), mult. (50->16), div. (0->0), fcn. (35->6), ass. (0->10)
	t104 = -pkin(3) + r_i_i_C(2);
	t103 = r_i_i_C(3) + qJ(4);
	t97 = qJ(2) + pkin(8);
	t96 = cos(t97);
	t102 = qJD(2) * t96;
	t95 = sin(t97);
	t101 = qJD(4) * t96 + (-cos(qJ(2)) * pkin(2) - t103 * t95 + t104 * t96) * qJD(2);
	t99 = cos(pkin(7));
	t98 = sin(pkin(7));
	t1 = [0, t101 * t99, 0, t99 * t102, 0; 0, t101 * t98, 0, t98 * t102, 0; 0, t95 * qJD(4) + (-sin(qJ(2)) * pkin(2) + t103 * t96 + t104 * t95) * qJD(2), 0, qJD(2) * t95, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:53
	% EndTime: 2019-10-24 10:22:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (82->23), mult. (142->45), div. (0->0), fcn. (107->8), ass. (0->21)
	t139 = sin(pkin(7));
	t141 = sin(qJ(5));
	t156 = t139 * t141;
	t142 = cos(qJ(5));
	t155 = t139 * t142;
	t140 = cos(pkin(7));
	t154 = t140 * t141;
	t153 = t140 * t142;
	t138 = qJ(2) + pkin(8);
	t136 = sin(t138);
	t152 = qJD(2) * t136;
	t137 = cos(t138);
	t151 = qJD(2) * t137;
	t150 = qJD(5) * t137;
	t149 = -pkin(3) - pkin(6) - r_i_i_C(3);
	t148 = r_i_i_C(1) * t142 - r_i_i_C(2) * t141;
	t147 = r_i_i_C(1) * t141 + r_i_i_C(2) * t142 + qJ(4);
	t146 = t148 * t151;
	t145 = t148 * qJD(5) + qJD(4);
	t144 = t145 * t137 + (-cos(qJ(2)) * pkin(2) - t147 * t136 + t149 * t137) * qJD(2);
	t1 = [0, t144 * t140, 0, t140 * t151, t140 * t146 + ((-t136 * t154 - t155) * r_i_i_C(1) + (-t136 * t153 + t156) * r_i_i_C(2)) * qJD(5); 0, t144 * t139, 0, t139 * t151, t139 * t146 + ((-t136 * t156 + t153) * r_i_i_C(1) + (-t136 * t155 - t154) * r_i_i_C(2)) * qJD(5); 0, t145 * t136 + (-sin(qJ(2)) * pkin(2) + t147 * t137 + t149 * t136) * qJD(2), 0, t152, (-t141 * t152 + t142 * t150) * r_i_i_C(2) + (t141 * t150 + t142 * t152) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end