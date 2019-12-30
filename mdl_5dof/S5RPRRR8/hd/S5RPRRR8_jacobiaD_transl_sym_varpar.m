% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
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
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (36->15), mult. (88->23), div. (0->0), fcn. (74->4), ass. (0->14)
	t87 = -pkin(1) - pkin(2);
	t78 = sin(qJ(1));
	t86 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t85 = qJD(1) * t80;
	t84 = qJD(3) * t78;
	t83 = qJD(3) * t80;
	t77 = sin(qJ(3));
	t79 = cos(qJ(3));
	t71 = -t77 * t84 - t79 * t83 + (t77 * t78 + t79 * t80) * qJD(1);
	t72 = (t84 - t86) * t79 + (-t83 + t85) * t77;
	t82 = -t72 * r_i_i_C(1) - t71 * r_i_i_C(2);
	t81 = t71 * r_i_i_C(1) - t72 * r_i_i_C(2);
	t1 = [t80 * qJD(2) + (-qJ(2) * t78 + t87 * t80) * qJD(1) - t81, t85, t81, 0, 0; t78 * qJD(2) + (qJ(2) * t80 + t87 * t78) * qJD(1) - t82, t86, t82, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:49:01
	% EndTime: 2019-12-29 17:49:01
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (123->23), mult. (308->39), div. (0->0), fcn. (276->6), ass. (0->21)
	t80 = qJD(1) - qJD(3);
	t68 = sin(qJ(3));
	t69 = sin(qJ(1));
	t70 = cos(qJ(3));
	t71 = cos(qJ(1));
	t47 = -t69 * t68 - t71 * t70;
	t45 = t80 * t47;
	t48 = t71 * t68 - t69 * t70;
	t46 = t80 * t48;
	t56 = sin(qJ(4));
	t57 = cos(qJ(4));
	t60 = r_i_i_C(1) * t57 - r_i_i_C(2) * t56 + pkin(3);
	t72 = pkin(7) + r_i_i_C(3);
	t77 = (r_i_i_C(1) * t56 + r_i_i_C(2) * t57) * qJD(4);
	t79 = -t72 * t45 - t60 * t46 - t47 * t77;
	t78 = -t60 * t45 + t72 * t46 + t48 * t77;
	t76 = qJ(2) * qJD(1);
	t73 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t67 = qJD(4) * t56;
	t66 = qJD(4) * t57;
	t1 = [-t69 * t76 + t73 * t71 - t78, qJD(1) * t71, t78, (-t46 * t57 - t47 * t67) * r_i_i_C(2) + (-t46 * t56 + t47 * t66) * r_i_i_C(1), 0; t73 * t69 + t71 * t76 - t79, qJD(1) * t69, t79, (t45 * t57 - t48 * t67) * r_i_i_C(2) + (t45 * t56 + t48 * t66) * r_i_i_C(1), 0; 0, 0, 0, t77, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:56
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (231->35), mult. (428->49), div. (0->0), fcn. (387->8), ass. (0->36)
	t103 = r_i_i_C(3) + pkin(8) + pkin(7);
	t110 = qJD(1) - qJD(3);
	t100 = sin(qJ(1));
	t101 = cos(qJ(3));
	t102 = cos(qJ(1));
	t99 = sin(qJ(3));
	t64 = -t100 * t99 - t102 * t101;
	t60 = t110 * t64;
	t65 = -t100 * t101 + t102 * t99;
	t61 = t110 * t65;
	t78 = qJD(4) + qJD(5);
	t79 = qJ(4) + qJ(5);
	t77 = cos(t79);
	t104 = r_i_i_C(2) * t77;
	t76 = sin(t79);
	t88 = r_i_i_C(1) * t76 + t104;
	t80 = sin(qJ(4));
	t95 = pkin(4) * qJD(4);
	t93 = t80 * t95;
	t83 = t88 * t78 + t93;
	t81 = cos(qJ(4));
	t87 = t81 * pkin(4) + r_i_i_C(1) * t77 - r_i_i_C(2) * t76 + pkin(3);
	t112 = -t103 * t60 - t87 * t61 - t83 * t64;
	t111 = t103 * t61 - t87 * t60 + t83 * t65;
	t109 = qJD(1) * t100;
	t108 = qJD(1) * t102;
	t105 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t98 = t76 * t78;
	t97 = t77 * t78;
	t96 = r_i_i_C(1) * t98 + r_i_i_C(2) * t97;
	t94 = r_i_i_C(2) * t98;
	t86 = r_i_i_C(1) * t97 + t81 * t95;
	t85 = -pkin(4) * t80 - t88;
	t63 = t64 * t94;
	t62 = t65 * t94;
	t1 = [-qJ(2) * t109 + t105 * t102 - t111, t108, t111, t85 * t61 + t86 * t64 - t63, -t61 * t104 - t63 + (-t61 * t76 + t64 * t97) * r_i_i_C(1); qJ(2) * t108 + t105 * t100 - t112, t109, t112, -t85 * t60 + t86 * t65 - t62, t60 * t104 - t62 + (t60 * t76 + t65 * t97) * r_i_i_C(1); 0, 0, 0, t93 + t96, t96;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end