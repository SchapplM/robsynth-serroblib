% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
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
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = r_i_i_C(3) + qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t34 * qJD(2) + (-t36 * t34 + t37 * t35) * qJD(1), qJD(1) * t34, 0, 0, 0; -t35 * qJD(2) + (t37 * t34 + t36 * t35) * qJD(1), -qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (21->14), mult. (62->21), div. (0->0), fcn. (48->6), ass. (0->11)
	t71 = sin(pkin(9));
	t72 = sin(pkin(8));
	t73 = cos(pkin(9));
	t83 = (r_i_i_C(3) + qJ(3)) * t72 + (r_i_i_C(1) * t73 - r_i_i_C(2) * t71 + pkin(2)) * cos(pkin(8)) + pkin(1);
	t75 = sin(qJ(1));
	t81 = qJD(1) * t75;
	t76 = cos(qJ(1));
	t80 = qJD(1) * t76;
	t79 = t72 * qJD(3);
	t77 = t71 * r_i_i_C(1) + t73 * r_i_i_C(2) + qJ(2);
	t1 = [0, 0, 0, 0, 0; t76 * t79 + t75 * qJD(2) + (-t83 * t75 + t77 * t76) * qJD(1), t81, t72 * t80, 0, 0; t75 * t79 - t76 * qJD(2) + (t77 * t75 + t83 * t76) * qJD(1), -t80, t72 * t81, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (84->29), mult. (140->49), div. (0->0), fcn. (118->8), ass. (0->23)
	t121 = sin(pkin(8));
	t122 = cos(pkin(8));
	t138 = (r_i_i_C(3) + pkin(6) + qJ(3)) * t121 + (cos(pkin(9)) * pkin(3) + pkin(2)) * t122 + pkin(1);
	t124 = sin(qJ(1));
	t136 = t122 * t124;
	t125 = cos(qJ(1));
	t135 = t122 * t125;
	t134 = qJD(1) * t124;
	t133 = qJD(1) * t125;
	t132 = qJD(3) * t121;
	t130 = sin(pkin(9)) * pkin(3) + qJ(2);
	t119 = pkin(9) + qJ(4);
	t117 = sin(t119);
	t118 = cos(t119);
	t129 = t117 * t124 + t118 * t135;
	t128 = t117 * t125 - t118 * t136;
	t127 = -t117 * t135 + t118 * t124;
	t126 = t117 * t136 + t118 * t125;
	t115 = t129 * qJD(1) - t126 * qJD(4);
	t114 = t127 * qJD(1) + t128 * qJD(4);
	t113 = t128 * qJD(1) + t127 * qJD(4);
	t112 = t126 * qJD(1) - t129 * qJD(4);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t118 + r_i_i_C(2) * t117) * t121 * qJD(4), 0; t125 * t132 + t113 * r_i_i_C(1) + t112 * r_i_i_C(2) + t124 * qJD(2) + (-t138 * t124 + t130 * t125) * qJD(1), t134, t121 * t133, t114 * r_i_i_C(1) - t115 * r_i_i_C(2), 0; t124 * t132 + t115 * r_i_i_C(1) + t114 * r_i_i_C(2) - t125 * qJD(2) + (t130 * t124 + t138 * t125) * qJD(1), -t133, t121 * t134, -t112 * r_i_i_C(1) + t113 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (233->46), mult. (239->69), div. (0->0), fcn. (198->10), ass. (0->34)
	t161 = pkin(9) + qJ(4);
	t158 = cos(t161);
	t163 = sin(pkin(8));
	t164 = cos(pkin(8));
	t185 = (r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3)) * t163 + (pkin(4) * t158 + cos(pkin(9)) * pkin(3) + pkin(2)) * t164 + pkin(1);
	t183 = pkin(4) * qJD(4);
	t165 = sin(qJ(1));
	t182 = t164 * t165;
	t166 = cos(qJ(1));
	t181 = t164 * t166;
	t157 = sin(t161);
	t180 = qJ(2) + pkin(4) * t157 + sin(pkin(9)) * pkin(3);
	t162 = qJD(4) + qJD(5);
	t159 = qJ(5) + t161;
	t155 = sin(t159);
	t156 = cos(t159);
	t169 = t155 * t182 + t156 * t166;
	t172 = t155 * t165 + t156 * t181;
	t148 = t169 * qJD(1) - t162 * t172;
	t170 = -t155 * t181 + t156 * t165;
	t171 = t155 * t166 - t156 * t182;
	t149 = qJD(1) * t171 + t162 * t170;
	t179 = -t148 * r_i_i_C(1) + t149 * r_i_i_C(2);
	t150 = qJD(1) * t170 + t162 * t171;
	t151 = qJD(1) * t172 - t169 * t162;
	t178 = t150 * r_i_i_C(1) - t151 * r_i_i_C(2);
	t177 = qJD(1) * t165;
	t176 = qJD(1) * t166;
	t175 = qJD(3) * t163;
	t174 = r_i_i_C(1) * t156 * t162;
	t168 = -t157 * t181 + t158 * t165;
	t167 = -t157 * t182 - t158 * t166;
	t152 = t163 * t162 * t155 * r_i_i_C(2);
	t1 = [0, 0, 0, t152 + (-t158 * t183 - t174) * t163, -t163 * t174 + t152; t166 * t175 + t149 * r_i_i_C(1) + t148 * r_i_i_C(2) + t165 * qJD(2) + t168 * t183 + (-t185 * t165 + t180 * t166) * qJD(1), t177, t163 * t176, ((t157 * t166 - t158 * t182) * qJD(4) + t168 * qJD(1)) * pkin(4) + t178, t178; t165 * t175 + t151 * r_i_i_C(1) + t150 * r_i_i_C(2) - t166 * qJD(2) + t167 * t183 + (t180 * t165 + t185 * t166) * qJD(1), -t176, t163 * t177, ((t157 * t165 + t158 * t181) * qJD(4) + t167 * qJD(1)) * pkin(4) + t179, t179;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end