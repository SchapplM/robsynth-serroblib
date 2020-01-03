% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t19 = pkin(1) * qJD(1);
	t16 = qJ(1) + qJ(2);
	t13 = sin(t16);
	t14 = cos(t16);
	t15 = qJD(1) + qJD(2);
	t18 = (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t15;
	t17 = (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t15;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t19 + t17, t17, 0, 0, 0; cos(qJ(1)) * t19 + t18, t18, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t25 = pkin(1) * qJD(1);
	t22 = qJ(1) + qJ(2);
	t19 = pkin(8) + t22;
	t17 = sin(t19);
	t18 = cos(t19);
	t21 = qJD(1) + qJD(2);
	t24 = (cos(t22) * pkin(2) + r_i_i_C(1) * t18 - r_i_i_C(2) * t17) * t21;
	t23 = (-pkin(2) * sin(t22) - r_i_i_C(1) * t17 - r_i_i_C(2) * t18) * t21;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t25 + t23, t23, 0, 0, 0; cos(qJ(1)) * t25 + t24, t24, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (109->16), mult. (66->20), div. (0->0), fcn. (40->8), ass. (0->15)
	t69 = qJD(1) + qJD(2);
	t83 = pkin(2) * t69;
	t82 = pkin(3) + r_i_i_C(1) * cos(pkin(9));
	t81 = -qJD(4) - r_i_i_C(2) * t69 * sin(pkin(9));
	t70 = qJ(1) + qJ(2);
	t67 = pkin(8) + t70;
	t65 = sin(t67);
	t79 = t69 * t65;
	t66 = cos(t67);
	t78 = t69 * t66;
	t77 = pkin(1) * qJD(1);
	t76 = qJ(4) * t69;
	t74 = t65 * t76 + r_i_i_C(3) * t79 + cos(t70) * t83 + t81 * t66 + t82 * t78;
	t73 = -sin(t70) * t83 + r_i_i_C(3) * t78 + t66 * t76 + (-t69 * t82 - t81) * t65;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t77 + t73, t73, 0, t79, 0; cos(qJ(1)) * t77 + t74, t74, 0, -t78, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (276->30), mult. (208->47), div. (0->0), fcn. (168->10), ass. (0->26)
	t166 = sin(pkin(9));
	t167 = cos(pkin(9));
	t185 = pkin(4) * t167 + (pkin(7) + r_i_i_C(3)) * t166 + pkin(3);
	t181 = pkin(1) * qJD(1);
	t165 = qJ(1) + qJ(2);
	t162 = pkin(8) + t165;
	t160 = sin(t162);
	t164 = qJD(1) + qJD(2);
	t180 = t164 * t160;
	t161 = cos(t162);
	t179 = t164 * t161;
	t168 = sin(qJ(5));
	t178 = t167 * t168;
	t169 = cos(qJ(5));
	t177 = t167 * t169;
	t175 = t160 * t168 + t161 * t177;
	t174 = t160 * t169 - t161 * t178;
	t173 = -t160 * t177 + t161 * t168;
	t172 = t160 * t178 + t161 * t169;
	t150 = qJD(5) * t173 + t164 * t174;
	t151 = -qJD(5) * t172 + t164 * t175;
	t171 = pkin(2) * t164 * cos(t165) + t151 * r_i_i_C(1) + t150 * r_i_i_C(2) + qJ(4) * t180 - t161 * qJD(4) + t185 * t179;
	t148 = -qJD(5) * t175 + t164 * t172;
	t149 = qJD(5) * t174 + t164 * t173;
	t170 = t148 * r_i_i_C(2) + t149 * r_i_i_C(1) + qJ(4) * t179 + t160 * qJD(4) + (-pkin(2) * sin(t165) - t185 * t160) * t164;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t169 + r_i_i_C(2) * t168) * t166 * qJD(5); -sin(qJ(1)) * t181 + t170, t170, 0, t180, t150 * r_i_i_C(1) - t151 * r_i_i_C(2); cos(qJ(1)) * t181 + t171, t171, 0, -t179, -t148 * r_i_i_C(1) + t149 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end