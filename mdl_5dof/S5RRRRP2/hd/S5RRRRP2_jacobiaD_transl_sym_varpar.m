% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
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
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
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
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (93->19), mult. (104->33), div. (0->0), fcn. (64->6), ass. (0->19)
	t61 = r_i_i_C(3) + pkin(7);
	t42 = qJ(1) + qJ(2);
	t39 = sin(t42);
	t40 = cos(t42);
	t44 = cos(qJ(3));
	t53 = qJD(3) * t44;
	t41 = qJD(1) + qJD(2);
	t43 = sin(qJ(3));
	t57 = t41 * t43;
	t60 = t39 * t57 - t40 * t53;
	t59 = t39 * t53 + t40 * t57;
	t56 = t41 * t44;
	t55 = pkin(1) * qJD(1);
	t54 = qJD(3) * t43;
	t50 = t39 * t54;
	t47 = t39 * t56 + t40 * t54;
	t46 = r_i_i_C(1) * t50 + ((-r_i_i_C(1) * t44 - pkin(2)) * t40 - t61 * t39) * t41 + t59 * r_i_i_C(2);
	t45 = -t47 * r_i_i_C(1) + t60 * r_i_i_C(2) + (-pkin(2) * t39 + t40 * t61) * t41;
	t1 = [-cos(qJ(1)) * t55 + t46, t46, t60 * r_i_i_C(1) + t47 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t55 + t45, t45, (-t40 * t56 + t50) * r_i_i_C(2) - t59 * r_i_i_C(1), 0, 0; 0, 0, (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (195->32), mult. (162->44), div. (0->0), fcn. (103->8), ass. (0->34)
	t62 = qJ(1) + qJ(2);
	t58 = cos(t62);
	t60 = qJD(1) + qJD(2);
	t83 = t58 * t60;
	t61 = qJ(3) + qJ(4);
	t55 = sin(t61);
	t56 = sin(t62);
	t57 = cos(t61);
	t59 = qJD(3) + qJD(4);
	t84 = t57 * t59;
	t89 = t55 * t83 + t56 * t84;
	t63 = sin(qJ(3));
	t88 = pkin(3) * t63;
	t87 = r_i_i_C(2) * t57;
	t86 = t55 * t59;
	t85 = t56 * t60;
	t82 = t60 * (-pkin(8) - pkin(7));
	t81 = pkin(1) * qJD(1);
	t64 = cos(qJ(3));
	t80 = qJD(3) * t64;
	t79 = r_i_i_C(1) * t84;
	t78 = t60 * t87;
	t77 = t56 * t86;
	t76 = t55 * t85;
	t73 = qJD(3) * t88;
	t72 = -pkin(3) * t64 - r_i_i_C(1) * t57 - pkin(2);
	t71 = -r_i_i_C(1) * t55 - t87;
	t70 = t71 * t59;
	t69 = r_i_i_C(1) * t76 + t56 * t78 + (r_i_i_C(2) * t86 - t79) * t58;
	t68 = t70 - t73;
	t67 = t72 * t83 + r_i_i_C(1) * t77 + (-r_i_i_C(3) * t60 + t73 + t82) * t56 + t89 * r_i_i_C(2);
	t66 = r_i_i_C(2) * t76 + r_i_i_C(3) * t83 + t72 * t85 + (t68 - t82) * t58;
	t43 = r_i_i_C(2) * t77;
	t1 = [-cos(qJ(1)) * t81 + t67, t67, (-t58 * t80 + t63 * t85) * pkin(3) + t69, t69, 0; -sin(qJ(1)) * t81 + t66, t66, t43 + (-pkin(3) * t80 - t79) * t56 + (t71 - t88) * t83, -t89 * r_i_i_C(1) - t58 * t78 + t43, 0; 0, 0, t68, t70, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (263->42), mult. (202->54), div. (0->0), fcn. (131->8), ass. (0->39)
	t98 = -pkin(4) - r_i_i_C(1);
	t71 = qJ(3) + qJ(4);
	t66 = cos(t71);
	t99 = t98 * t66;
	t64 = sin(t71);
	t97 = pkin(4) * t64;
	t96 = r_i_i_C(2) * t66;
	t72 = qJ(1) + qJ(2);
	t65 = sin(t72);
	t69 = qJD(3) + qJD(4);
	t95 = t65 * t69;
	t67 = cos(t72);
	t93 = t67 * t69;
	t70 = qJD(1) + qJD(2);
	t92 = t70 * t65;
	t91 = t70 * t67;
	t90 = pkin(1) * qJD(1);
	t89 = pkin(3) * qJD(3);
	t88 = t65 * t96;
	t87 = r_i_i_C(2) * t91;
	t86 = t64 * t95;
	t85 = t64 * t92;
	t84 = t66 * t93;
	t83 = t64 * r_i_i_C(2) * t93 + r_i_i_C(1) * t85 + t70 * t88;
	t73 = sin(qJ(3));
	t82 = t73 * t89;
	t81 = t98 * t64;
	t74 = cos(qJ(3));
	t80 = -t74 * pkin(3) - pkin(2) + t99;
	t79 = t69 * t99 - t74 * t89;
	t78 = -r_i_i_C(1) * t64 - t96;
	t77 = (t81 - t96) * t69;
	t49 = -t69 * t97 - t82;
	t68 = qJ(5) + pkin(7) + pkin(8);
	t76 = -t65 * t49 + r_i_i_C(1) * t86 + t69 * t88 + t64 * t87 + qJD(5) * t67 + (t80 * t67 + (-r_i_i_C(3) - t68) * t65) * t70;
	t75 = r_i_i_C(2) * t85 + r_i_i_C(3) * t91 + t65 * qJD(5) + t80 * t92 + (t68 * t70 + t78 * t69 + t49) * t67;
	t60 = -t73 * pkin(3) - t97;
	t51 = r_i_i_C(2) * t86;
	t1 = [-cos(qJ(1)) * t90 + t76, t76, -t60 * t92 + t79 * t67 + t83, -r_i_i_C(1) * t84 + (-t84 + t85) * pkin(4) + t83, t91; -sin(qJ(1)) * t90 + t75, t75, t51 + t79 * t65 + (t60 + t78) * t91, t51 + t81 * t91 + (t98 * t95 - t87) * t66, t92; 0, 0, t77 - t82, t77, 0;];
	JaD_transl = t1;
end