% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR4
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t68 = sin(qJ(3));
	t69 = cos(qJ(3));
	t70 = qJD(3) * (r_i_i_C(1) * t68 + r_i_i_C(2) * t69);
	t77 = -pkin(6) - r_i_i_C(3);
	t76 = qJD(1) * t68;
	t75 = qJD(1) * t69;
	t74 = qJD(3) * t68;
	t73 = qJD(3) * t69;
	t71 = r_i_i_C(1) * t69 - r_i_i_C(2) * t68 + pkin(2);
	t67 = qJ(1) + pkin(8);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [0, 0, -t70, 0, 0; t66 * t70 + (sin(qJ(1)) * pkin(1) + t77 * t66 + t71 * t65) * qJD(1), 0, (-t65 * t74 + t66 * t75) * r_i_i_C(2) + (t65 * t73 + t66 * t76) * r_i_i_C(1), 0, 0; t65 * t70 + (-cos(qJ(1)) * pkin(1) + t77 * t65 - t71 * t66) * qJD(1), 0, (t65 * t75 + t66 * t74) * r_i_i_C(2) + (t65 * t76 - t66 * t73) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (85->20), mult. (94->27), div. (0->0), fcn. (61->8), ass. (0->14)
	t79 = qJ(3) + pkin(9);
	t75 = sin(t79);
	t77 = cos(t79);
	t86 = r_i_i_C(1) * t77 - r_i_i_C(2) * t75 + cos(qJ(3)) * pkin(3);
	t92 = qJD(3) * t86;
	t87 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t75 + r_i_i_C(2) * t77;
	t84 = qJD(3) * t87;
	t90 = -r_i_i_C(3) - qJ(4) - pkin(6);
	t88 = pkin(2) + t86;
	t85 = qJD(1) * t87;
	t80 = qJ(1) + pkin(8);
	t78 = cos(t80);
	t76 = sin(t80);
	t1 = [0, 0, -t84, 0, 0; -t76 * qJD(4) + t78 * t84 + (sin(qJ(1)) * pkin(1) + t90 * t78 + t88 * t76) * qJD(1), 0, t76 * t92 + t78 * t85, -qJD(1) * t76, 0; t78 * qJD(4) + t76 * t84 + (-cos(qJ(1)) * pkin(1) + t90 * t76 - t88 * t78) * qJD(1), 0, t76 * t85 - t78 * t92, qJD(1) * t78, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:54:47
	% EndTime: 2019-12-05 17:54:47
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (174->30), mult. (136->38), div. (0->0), fcn. (88->10), ass. (0->26)
	t112 = qJ(3) + pkin(9);
	t109 = qJ(5) + t112;
	t103 = sin(t109);
	t104 = cos(t109);
	t136 = r_i_i_C(1) * t103 + r_i_i_C(2) * t104;
	t135 = t136 * qJD(1);
	t102 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t112);
	t111 = qJD(3) + qJD(5);
	t117 = t136 * t111;
	t134 = t102 * qJD(3) - t117;
	t132 = r_i_i_C(1) * t104;
	t131 = r_i_i_C(2) * t103;
	t129 = -r_i_i_C(3) - pkin(7) - qJ(4) - pkin(6);
	t113 = qJ(1) + pkin(8);
	t106 = sin(t113);
	t128 = qJD(1) * t106;
	t108 = cos(t113);
	t127 = qJD(1) * t108;
	t124 = t111 * t131;
	t126 = t135 * t106 + t108 * t124;
	t125 = t111 * t132;
	t123 = t106 * t125 + t135 * t108;
	t120 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t112);
	t118 = pkin(2) - t120 - t131 + t132;
	t96 = t120 * qJD(3);
	t1 = [0, 0, t134, 0, -t117; -t106 * qJD(4) - t134 * t108 + (sin(qJ(1)) * pkin(1) + t129 * t108 + t118 * t106) * qJD(1), 0, -t102 * t127 + (-t96 - t124) * t106 + t123, -t128, -t106 * t124 + t123; t108 * qJD(4) - t134 * t106 + (-cos(qJ(1)) * pkin(1) + t129 * t106 - t118 * t108) * qJD(1), 0, -t102 * t128 + (t96 - t125) * t108 + t126, t127, -t108 * t125 + t126;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end