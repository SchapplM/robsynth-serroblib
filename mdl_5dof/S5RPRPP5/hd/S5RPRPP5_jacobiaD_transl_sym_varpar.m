% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(6) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t22 * t18 + t24 * t20) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0; t21 * t18 + (t24 * t18 + t22 * t20) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0; 0, 0, -t29, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:36
	% EndTime: 2019-12-29 16:43:36
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (48->21), mult. (142->35), div. (0->0), fcn. (98->4), ass. (0->17)
	t134 = sin(qJ(3));
	t136 = cos(qJ(3));
	t146 = r_i_i_C(3) + qJ(4);
	t147 = pkin(3) + r_i_i_C(1);
	t151 = -(t146 * t134 + t147 * t136) * qJD(3) + t136 * qJD(4);
	t150 = t147 * qJD(3) - qJD(4);
	t148 = t147 * t134 - t146 * t136 + qJ(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(3) * t135;
	t142 = qJD(3) * t137;
	t140 = qJD(1) * t147;
	t139 = qJD(1) * t146;
	t138 = qJD(2) + (-pkin(1) - pkin(6) - r_i_i_C(2)) * qJD(1) - t151;
	t1 = [t138 * t137 - t148 * t145, t144, (t137 * t140 + t146 * t143) * t136 + (-t150 * t135 + t137 * t139) * t134, t134 * t143 - t136 * t144, 0; t138 * t135 + t148 * t144, t145, (t135 * t140 - t146 * t142) * t136 + (t135 * t139 + t150 * t137) * t134, -t134 * t142 - t136 * t145, 0; 0, 0, t151, qJD(3) * t136, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:35
	% EndTime: 2019-12-29 16:43:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (62->24), mult. (176->36), div. (0->0), fcn. (121->4), ass. (0->17)
	t20 = sin(qJ(3));
	t22 = cos(qJ(3));
	t27 = pkin(3) + pkin(4) + r_i_i_C(1);
	t32 = r_i_i_C(2) + qJ(4);
	t37 = -(t32 * t20 + t27 * t22) * qJD(3) + t22 * qJD(4);
	t36 = -qJD(5) + (t27 * t20 - t32 * t22 + qJ(2)) * qJD(1);
	t35 = t27 * qJD(3) - qJD(4);
	t21 = sin(qJ(1));
	t19 = qJD(1) * t21;
	t23 = cos(qJ(1));
	t31 = qJD(1) * t23;
	t30 = qJD(3) * t21;
	t29 = qJD(3) * t23;
	t26 = qJD(1) * t32;
	t25 = qJD(1) * t27;
	t24 = qJD(2) + (-pkin(1) - pkin(6) + r_i_i_C(3) + qJ(5)) * qJD(1) - t37;
	t1 = [-t36 * t21 + t24 * t23, t31, (t23 * t25 + t32 * t30) * t22 + (-t35 * t21 + t23 * t26) * t20, t20 * t30 - t22 * t31, t19; t24 * t21 + t36 * t23, t19, (t21 * t25 - t32 * t29) * t22 + (t21 * t26 + t35 * t23) * t20, -t22 * t19 - t20 * t29, -t31; 0, 0, t37, qJD(3) * t22, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end