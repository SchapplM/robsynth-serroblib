% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR3
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
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
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
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
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->11), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t42 = r_i_i_C(3) + qJ(3);
	t41 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(2);
	t38 = qJ(1) + pkin(8);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [0, 0, 0, 0, 0; t36 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t42 * t37 - t41 * t36) * qJD(1), 0, qJD(1) * t36, 0, 0; -t37 * qJD(3) + (cos(qJ(1)) * pkin(1) + t42 * t36 + t41 * t37) * qJD(1), 0, -qJD(1) * t37, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:14
	% EndTime: 2020-01-03 11:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (72->22), mult. (74->34), div. (0->0), fcn. (48->7), ass. (0->14)
	t85 = r_i_i_C(3) + pkin(6) + qJ(3);
	t77 = qJ(1) + pkin(8);
	t73 = sin(t77);
	t84 = qJD(1) * t73;
	t75 = cos(t77);
	t83 = qJD(1) * t75;
	t82 = qJD(4) * t73;
	t81 = qJD(4) * t75;
	t76 = pkin(9) + qJ(4);
	t72 = sin(t76);
	t74 = cos(t76);
	t80 = r_i_i_C(1) * t74 - r_i_i_C(2) * t72 + cos(pkin(9)) * pkin(3) + pkin(2);
	t79 = (-r_i_i_C(1) * t72 - r_i_i_C(2) * t74) * qJD(4);
	t1 = [0, 0, 0, t79, 0; t73 * qJD(3) + t75 * t79 + (-sin(qJ(1)) * pkin(1) + t85 * t75 - t80 * t73) * qJD(1), 0, t84, (t72 * t82 - t74 * t83) * r_i_i_C(2) + (-t72 * t83 - t74 * t82) * r_i_i_C(1), 0; -t75 * qJD(3) + t73 * t79 + (cos(qJ(1)) * pkin(1) + t85 * t73 + t80 * t75) * qJD(1), 0, -t83, (-t72 * t81 - t74 * t84) * r_i_i_C(2) + (-t72 * t84 + t74 * t81) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:29:15
	% EndTime: 2020-01-03 11:29:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (167->33), mult. (124->44), div. (0->0), fcn. (81->9), ass. (0->29)
	t103 = pkin(9) + qJ(4);
	t101 = qJ(5) + t103;
	t95 = sin(t101);
	t121 = r_i_i_C(1) * t95;
	t96 = cos(t101);
	t120 = r_i_i_C(2) * t96;
	t119 = r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3);
	t118 = pkin(4) * qJD(4);
	t104 = qJD(4) + qJD(5);
	t117 = t104 * t95;
	t116 = t104 * t96;
	t105 = qJ(1) + pkin(8);
	t98 = sin(t105);
	t115 = qJD(1) * t98;
	t100 = cos(t105);
	t114 = qJD(1) * t100;
	t113 = r_i_i_C(1) * t116;
	t112 = r_i_i_C(2) * t117;
	t99 = cos(t103);
	t111 = t99 * t118;
	t110 = -t120 - t121;
	t109 = r_i_i_C(1) * t96 - r_i_i_C(2) * t95 + pkin(4) * t99 + cos(pkin(9)) * pkin(3) + pkin(2);
	t108 = t110 * t104;
	t97 = sin(t103);
	t107 = qJD(1) * (-pkin(4) * t97 + t110);
	t106 = -t97 * t118 + t108;
	t93 = t100 * t113;
	t92 = t98 * t112;
	t1 = [0, 0, 0, t106, t108; t98 * qJD(3) + t106 * t100 + (-sin(qJ(1)) * pkin(1) + t119 * t100 - t109 * t98) * qJD(1), 0, t115, t92 + (-t111 - t113) * t98 + t100 * t107, -t114 * t120 + t92 + (-t95 * t114 - t98 * t116) * r_i_i_C(1); -t100 * qJD(3) + t106 * t98 + (cos(qJ(1)) * pkin(1) + t119 * t98 + t109 * t100) * qJD(1), 0, -t114, t93 + (t111 - t112) * t100 + t98 * t107, -t115 * t121 + t93 + (-t100 * t117 - t96 * t115) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end