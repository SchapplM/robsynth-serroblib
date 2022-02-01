% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(1) + qJD(2);
	t55 = pkin(2) * t48;
	t54 = pkin(1) * qJD(1);
	t49 = qJ(1) + qJ(2);
	t47 = qJ(3) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t44 = qJD(3) + t48;
	t53 = (-r_i_i_C(1) * t43 + r_i_i_C(2) * t42) * t44;
	t52 = (-r_i_i_C(1) * t42 - r_i_i_C(2) * t43) * t44;
	t51 = -cos(t49) * t55 + t53;
	t50 = -sin(t49) * t55 + t52;
	t1 = [-cos(qJ(1)) * t54 + t51, t51, t53, 0, 0; -sin(qJ(1)) * t54 + t50, t50, t52, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:46
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (225->23), mult. (148->36), div. (0->0), fcn. (90->8), ass. (0->24)
	t73 = r_i_i_C(3) + pkin(8);
	t51 = qJ(1) + qJ(2);
	t49 = qJ(3) + t51;
	t44 = sin(t49);
	t45 = cos(t49);
	t53 = cos(qJ(4));
	t64 = qJD(4) * t53;
	t50 = qJD(1) + qJD(2);
	t46 = qJD(3) + t50;
	t52 = sin(qJ(4));
	t68 = t46 * t52;
	t72 = t44 * t68 - t45 * t64;
	t71 = t44 * t64 + t45 * t68;
	t70 = pkin(2) * t50;
	t67 = t46 * t53;
	t66 = pkin(1) * qJD(1);
	t65 = qJD(4) * t52;
	t61 = t44 * t65;
	t58 = t44 * t67 + t45 * t65;
	t57 = r_i_i_C(1) * t61 + ((-r_i_i_C(1) * t53 - pkin(3)) * t45 - t73 * t44) * t46 + t71 * r_i_i_C(2);
	t56 = -cos(t51) * t70 + t57;
	t55 = -t58 * r_i_i_C(1) + t72 * r_i_i_C(2) + (-pkin(3) * t44 + t45 * t73) * t46;
	t54 = -sin(t51) * t70 + t55;
	t1 = [-cos(qJ(1)) * t66 + t56, t56, t57, t72 * r_i_i_C(1) + t58 * r_i_i_C(2), 0; -sin(qJ(1)) * t66 + t54, t54, t55, (-t45 * t67 + t61) * r_i_i_C(2) - t71 * r_i_i_C(1), 0; 0, 0, 0, (-r_i_i_C(1) * t52 - r_i_i_C(2) * t53) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:46
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (371->36), mult. (214->48), div. (0->0), fcn. (135->10), ass. (0->40)
	t71 = qJ(1) + qJ(2);
	t67 = qJ(3) + t71;
	t60 = cos(t67);
	t69 = qJD(1) + qJD(2);
	t62 = qJD(3) + t69;
	t97 = t60 * t62;
	t59 = sin(t67);
	t70 = qJ(4) + qJ(5);
	t65 = cos(t70);
	t68 = qJD(4) + qJD(5);
	t93 = t65 * t68;
	t63 = sin(t70);
	t96 = t62 * t63;
	t102 = t59 * t93 + t60 * t96;
	t101 = pkin(2) * t69;
	t72 = sin(qJ(4));
	t100 = pkin(4) * t72;
	t99 = r_i_i_C(2) * t65;
	t98 = t59 * t62;
	t95 = t62 * (-pkin(9) - pkin(8));
	t94 = t63 * t68;
	t92 = pkin(1) * qJD(1);
	t73 = cos(qJ(4));
	t91 = qJD(4) * t73;
	t90 = r_i_i_C(1) * t93;
	t89 = t62 * t99;
	t88 = t59 * t96;
	t87 = t59 * t94;
	t84 = qJD(4) * t100;
	t83 = -t73 * pkin(4) - r_i_i_C(1) * t65 - pkin(3);
	t82 = -r_i_i_C(1) * t63 - t99;
	t81 = t82 * t68;
	t80 = r_i_i_C(1) * t88 + t59 * t89 + (r_i_i_C(2) * t94 - t90) * t60;
	t79 = t81 - t84;
	t78 = t83 * t97 + r_i_i_C(1) * t87 + (-r_i_i_C(3) * t62 + t84 + t95) * t59 + t102 * r_i_i_C(2);
	t77 = -cos(t71) * t101 + t78;
	t76 = r_i_i_C(2) * t88 + r_i_i_C(3) * t97 + t83 * t98 + (t79 - t95) * t60;
	t75 = -sin(t71) * t101 + t76;
	t52 = r_i_i_C(2) * t87;
	t1 = [-cos(qJ(1)) * t92 + t77, t77, t78, (-t60 * t91 + t72 * t98) * pkin(4) + t80, t80; -sin(qJ(1)) * t92 + t75, t75, t76, t52 + (-pkin(4) * t91 - t90) * t59 + (t82 - t100) * t97, -t102 * r_i_i_C(1) - t60 * t89 + t52; 0, 0, 0, t79, t81;];
	JaD_transl = t1;
end