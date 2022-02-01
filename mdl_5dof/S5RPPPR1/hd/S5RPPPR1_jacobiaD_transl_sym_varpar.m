% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(7);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t15 = qJ(1) + pkin(7);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (49->17), mult. (66->24), div. (0->0), fcn. (50->8), ass. (0->11)
	t115 = sin(pkin(8));
	t121 = qJD(1) * t115;
	t120 = t115 * qJD(4);
	t114 = sin(pkin(9));
	t116 = cos(pkin(9));
	t119 = t114 * r_i_i_C(1) + t116 * r_i_i_C(2) + qJ(3);
	t118 = -pkin(2) + (-r_i_i_C(3) - qJ(4)) * t115 + (-r_i_i_C(1) * t116 + r_i_i_C(2) * t114 - pkin(3)) * cos(pkin(8));
	t113 = qJ(1) + pkin(7);
	t112 = cos(t113);
	t111 = sin(t113);
	t1 = [-t111 * t120 + t112 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t119 * t111 + t118 * t112) * qJD(1), 0, qJD(1) * t112, -t111 * t121, 0; t112 * t120 + t111 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t119 * t112 + t118 * t111) * qJD(1), 0, qJD(1) * t111, t112 * t121, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (138->32), mult. (144->52), div. (0->0), fcn. (120->10), ass. (0->23)
	t154 = qJ(1) + pkin(7);
	t150 = sin(t154);
	t157 = cos(pkin(8));
	t168 = t150 * t157;
	t152 = cos(t154);
	t167 = t152 * t157;
	t156 = sin(pkin(8));
	t166 = qJD(1) * t156;
	t165 = qJD(4) * t156;
	t164 = pkin(4) * sin(pkin(9)) + qJ(3);
	t153 = pkin(9) + qJ(5);
	t149 = sin(t153);
	t151 = cos(t153);
	t163 = -t149 * t150 - t151 * t167;
	t162 = -t149 * t152 + t151 * t168;
	t161 = t149 * t167 - t150 * t151;
	t160 = t149 * t168 + t151 * t152;
	t159 = -(cos(pkin(9)) * pkin(4) + pkin(3)) * t157 - pkin(2) + (-r_i_i_C(3) - pkin(6) - qJ(4)) * t156;
	t147 = t163 * qJD(1) + t160 * qJD(5);
	t146 = t161 * qJD(1) + t162 * qJD(5);
	t145 = t162 * qJD(1) + t161 * qJD(5);
	t144 = t160 * qJD(1) + t163 * qJD(5);
	t1 = [-t150 * t165 + t147 * r_i_i_C(1) + t146 * r_i_i_C(2) + t152 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t164 * t150 + t159 * t152) * qJD(1), 0, qJD(1) * t152, -t150 * t166, t144 * r_i_i_C(1) + t145 * r_i_i_C(2); t152 * t165 - t145 * r_i_i_C(1) + t144 * r_i_i_C(2) + t150 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t164 * t152 + t159 * t150) * qJD(1), 0, qJD(1) * t150, t152 * t166, -t146 * r_i_i_C(1) + t147 * r_i_i_C(2); 0, 0, 0, 0, (-r_i_i_C(1) * t151 + r_i_i_C(2) * t149) * t156 * qJD(5);];
	JaD_transl = t1;
end