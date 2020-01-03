% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
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
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(8);
	t16 = qJ(3) + t18;
	t14 = sin(t16);
	t15 = cos(t16);
	t17 = qJD(1) + qJD(3);
	t20 = (r_i_i_C(1) * t15 - r_i_i_C(2) * t14) * t17;
	t19 = (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15) * t17;
	t1 = [0, 0, 0, 0, 0; t19 + (-sin(t18) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t19, 0, 0; (cos(t18) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:01
	% EndTime: 2020-01-03 11:37:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (101->16), mult. (62->20), div. (0->0), fcn. (38->8), ass. (0->13)
	t77 = pkin(3) + r_i_i_C(1) * cos(pkin(9));
	t65 = qJD(1) + qJD(3);
	t76 = qJD(4) + r_i_i_C(2) * t65 * sin(pkin(9));
	t66 = qJ(1) + pkin(8);
	t64 = qJ(3) + t66;
	t62 = sin(t64);
	t74 = t65 * t62;
	t63 = cos(t64);
	t73 = t65 * t63;
	t72 = qJ(4) * t65;
	t70 = r_i_i_C(3) * t73 + t62 * t76 + t63 * t72 - t74 * t77;
	t69 = r_i_i_C(3) * t74 + t62 * t72 - t63 * t76 + t73 * t77;
	t1 = [0, 0, 0, 0, 0; (-sin(t66) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t70, 0, t70, t74, 0; (cos(t66) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t69, 0, t69, -t73, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:37:02
	% EndTime: 2020-01-03 11:37:02
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (268->30), mult. (204->46), div. (0->0), fcn. (166->10), ass. (0->25)
	t162 = sin(pkin(9));
	t163 = cos(pkin(9));
	t180 = pkin(4) * t163 + (pkin(7) + r_i_i_C(3)) * t162 + pkin(3);
	t161 = qJ(1) + pkin(8);
	t159 = qJ(3) + t161;
	t157 = sin(t159);
	t160 = qJD(1) + qJD(3);
	t176 = t160 * t157;
	t158 = cos(t159);
	t175 = t160 * t158;
	t164 = sin(qJ(5));
	t174 = t163 * t164;
	t165 = cos(qJ(5));
	t173 = t163 * t165;
	t171 = t157 * t164 + t158 * t173;
	t170 = t157 * t165 - t158 * t174;
	t169 = -t157 * t173 + t158 * t164;
	t168 = t157 * t174 + t158 * t165;
	t148 = t169 * qJD(5) + t170 * t160;
	t149 = -t168 * qJD(5) + t171 * t160;
	t167 = t149 * r_i_i_C(1) + t148 * r_i_i_C(2) + qJ(4) * t176 - t158 * qJD(4) + t180 * t175;
	t146 = -t171 * qJD(5) + t168 * t160;
	t147 = t170 * qJD(5) + t169 * t160;
	t166 = t147 * r_i_i_C(1) + t146 * r_i_i_C(2) + qJ(4) * t175 + t157 * qJD(4) - t180 * t176;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t165 + r_i_i_C(2) * t164) * t162 * qJD(5); (-sin(t161) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t166, 0, t166, t176, t148 * r_i_i_C(1) - t149 * r_i_i_C(2); (cos(t161) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t167, 0, t167, -t175, -t146 * r_i_i_C(1) + t147 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end