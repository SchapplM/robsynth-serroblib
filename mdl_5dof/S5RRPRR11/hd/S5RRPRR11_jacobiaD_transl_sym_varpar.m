% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR11_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR11_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR11_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:16:27
	% EndTime: 2019-12-29 19:16:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:16:27
	% EndTime: 2019-12-29 19:16:27
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
	% StartTime: 2019-12-29 19:16:15
	% EndTime: 2019-12-29 19:16:15
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:16:29
	% EndTime: 2019-12-29 19:16:29
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t134 = sin(qJ(2));
	t136 = cos(qJ(2));
	t146 = r_i_i_C(3) + qJ(3);
	t148 = pkin(2) + r_i_i_C(1);
	t149 = t148 * t134 - t146 * t136;
	t150 = qJD(2) * t149 - t134 * qJD(3);
	t147 = pkin(6) + r_i_i_C(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(2) * t137;
	t141 = -t146 * t134 - t148 * t136;
	t139 = -pkin(1) + t141;
	t1 = [t150 * t135 + (-t147 * t135 + t139 * t137) * qJD(1), (-t146 * t143 + t148 * t145) * t134 + (-t146 * t145 + (-t148 * qJD(2) + qJD(3)) * t137) * t136, -t134 * t145 + t136 * t143, 0, 0; -t150 * t137 + (t139 * t135 + t147 * t137) * qJD(1), -t149 * t144 + (t141 * qJD(2) + qJD(3) * t136) * t135, t135 * qJD(2) * t136 + t134 * t144, 0, 0; 0, -t150, qJD(2) * t134, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:16:15
	% EndTime: 2019-12-29 19:16:16
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (122->38), mult. (376->64), div. (0->0), fcn. (325->6), ass. (0->33)
	t59 = sin(qJ(4));
	t60 = sin(qJ(2));
	t62 = cos(qJ(4));
	t63 = cos(qJ(2));
	t72 = t59 * t63 - t60 * t62;
	t90 = pkin(2) + pkin(3);
	t91 = -qJ(3) * t63 + t90 * t60;
	t94 = t91 * qJD(2) - t60 * qJD(3);
	t71 = t59 * t60 + t62 * t63;
	t93 = qJD(2) - qJD(4);
	t54 = t93 * t71;
	t84 = qJD(2) * t60;
	t92 = t72 * qJD(4) + t62 * t84;
	t61 = sin(qJ(1));
	t86 = qJD(1) * t61;
	t64 = cos(qJ(1));
	t85 = qJD(1) * t64;
	t83 = qJD(2) * t63;
	t82 = qJD(2) * t64;
	t80 = pkin(6) - pkin(7) - r_i_i_C(3);
	t76 = t63 * t82;
	t68 = qJD(1) * t72;
	t49 = t54 * t64 + t61 * t68;
	t67 = qJD(1) * t71;
	t50 = -t59 * t76 + t61 * t67 + t92 * t64;
	t75 = t49 * r_i_i_C(1) + t50 * r_i_i_C(2);
	t51 = -t61 * t54 + t64 * t68;
	t52 = t64 * t67 + (t59 * t83 - t92) * t61;
	t74 = -t51 * r_i_i_C(1) - t52 * r_i_i_C(2);
	t73 = -t93 * t72 * r_i_i_C(1) - t54 * r_i_i_C(2);
	t70 = -qJ(3) * t60 - t90 * t63;
	t66 = -pkin(1) + t70;
	t1 = [-t52 * r_i_i_C(1) + t51 * r_i_i_C(2) + t94 * t61 + (-t80 * t61 + t66 * t64) * qJD(1), (-qJ(3) * t82 + t90 * t86) * t60 + (-qJ(3) * t86 + (-t90 * qJD(2) + qJD(3)) * t64) * t63 - t75, -t60 * t86 + t76, t75, 0; -t50 * r_i_i_C(1) + t49 * r_i_i_C(2) - t94 * t64 + (t66 * t61 + t80 * t64) * qJD(1), -t91 * t85 + (qJD(2) * t70 + qJD(3) * t63) * t61 - t74, t60 * t85 + t61 * t83, t74, 0; 0, -t94 - t73, t84, t73, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:16:21
	% EndTime: 2019-12-29 19:16:21
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (361->53), mult. (612->76), div. (0->0), fcn. (520->8), ass. (0->44)
	t85 = qJ(4) + qJ(5);
	t82 = sin(t85);
	t83 = cos(t85);
	t87 = sin(qJ(2));
	t90 = cos(qJ(2));
	t106 = t82 * t90 - t83 * t87;
	t86 = sin(qJ(4));
	t107 = pkin(4) * t86 + qJ(3);
	t89 = cos(qJ(4));
	t124 = t89 * pkin(4) + pkin(2) + pkin(3);
	t97 = t107 * t90 - t124 * t87;
	t128 = t97 * qJD(2) + t87 * qJD(3);
	t105 = t82 * t87 + t83 * t90;
	t84 = qJD(4) + qJD(5);
	t127 = qJD(2) - t84;
	t76 = t127 * t105;
	t102 = qJD(1) * t106;
	t88 = sin(qJ(1));
	t91 = cos(qJ(1));
	t69 = t88 * t102 + t76 * t91;
	t101 = qJD(1) * t105;
	t114 = qJD(2) * t90;
	t108 = t91 * t114;
	t115 = qJD(2) * t87;
	t125 = t106 * t84 + t83 * t115;
	t70 = t88 * t101 - t82 * t108 + t125 * t91;
	t121 = t69 * r_i_i_C(1) + t70 * r_i_i_C(2);
	t71 = t91 * t102 - t88 * t76;
	t72 = t91 * t101 + (t82 * t114 - t125) * t88;
	t120 = -t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
	t119 = -t127 * t106 * r_i_i_C(1) - t76 * r_i_i_C(2);
	t126 = (pkin(6) - r_i_i_C(3) - pkin(8) - pkin(7)) * qJD(1);
	t118 = pkin(4) * qJD(4);
	t117 = qJD(1) * t88;
	t116 = qJD(1) * t91;
	t104 = t86 * t90 - t87 * t89;
	t103 = t86 * t87 + t89 * t90;
	t100 = t104 * qJD(4);
	t99 = -t107 * t87 - t124 * t90;
	t96 = qJD(1) * (-pkin(1) + t99);
	t95 = (qJD(2) - qJD(4)) * t103;
	t94 = t99 * qJD(2) + qJD(3) * t90 + t103 * t118;
	t93 = -t104 * t118 + t128;
	t1 = [-t72 * r_i_i_C(1) + t71 * r_i_i_C(2) + t91 * t96 + (pkin(4) * t100 - t126 - t128) * t88, -t117 * t97 + t94 * t91 - t121, -t87 * t117 + t108, (t104 * t117 + t95 * t91) * pkin(4) + t121, t121; -t70 * r_i_i_C(1) + t69 * r_i_i_C(2) + t88 * t96 + (t93 + t126) * t91, t97 * t116 + t94 * t88 - t120, t88 * t114 + t87 * t116, (-t104 * t116 + t95 * t88) * pkin(4) + t120, t120; 0, t93 - t119, t115, (-t104 * qJD(2) + t100) * pkin(4) + t119, t119;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end