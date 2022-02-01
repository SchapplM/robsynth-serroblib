% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
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
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
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
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
	% DurationCPUTime: 0.09s
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
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (158->22), mult. (140->27), div. (0->0), fcn. (91->8), ass. (0->21)
	t52 = qJD(1) + qJD(2);
	t53 = qJ(3) + pkin(9);
	t48 = sin(t53);
	t49 = cos(t53);
	t62 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49 + sin(qJ(3)) * pkin(3);
	t60 = t62 * qJD(3);
	t78 = -t60 - t52 * (-qJ(4) - pkin(7));
	t71 = r_i_i_C(2) * t48;
	t77 = t52 * t71 + qJD(4);
	t76 = -r_i_i_C(1) * t49 - cos(qJ(3)) * pkin(3);
	t54 = qJ(1) + qJ(2);
	t50 = sin(t54);
	t69 = t52 * t50;
	t51 = cos(t54);
	t68 = t52 * t51;
	t66 = pkin(1) * qJD(1);
	t63 = -pkin(2) + t76;
	t61 = qJD(3) * (t71 + t76);
	t59 = (t63 * t52 + t77) * t51 + (-r_i_i_C(3) * t52 - t78) * t50;
	t58 = r_i_i_C(3) * t68 + t77 * t50 + t78 * t51 + t63 * t69;
	t1 = [-cos(qJ(1)) * t66 + t59, t59, t51 * t61 + t62 * t69, t68, 0; -sin(qJ(1)) * t66 + t58, t58, t50 * t61 - t62 * t68, t69, 0; 0, 0, -t60, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:44:01
	% EndTime: 2022-01-20 11:44:01
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (273->36), mult. (188->48), div. (0->0), fcn. (122->10), ass. (0->34)
	t73 = qJ(3) + pkin(9);
	t59 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t73);
	t56 = t59 * qJD(3);
	t71 = qJD(3) + qJD(5);
	t67 = qJ(5) + t73;
	t63 = sin(t67);
	t64 = cos(t67);
	t95 = r_i_i_C(2) * t64;
	t80 = -r_i_i_C(1) * t63 - t95;
	t79 = t80 * t71;
	t98 = t79 + t56;
	t74 = qJ(1) + qJ(2);
	t69 = cos(t74);
	t72 = qJD(1) + qJD(2);
	t92 = t72 * t69;
	t68 = sin(t74);
	t94 = t68 * t71;
	t97 = t63 * t92 + t64 * t94;
	t96 = r_i_i_C(1) * t64;
	t93 = t72 * t68;
	t91 = pkin(1) * qJD(1);
	t90 = t71 * t96;
	t89 = t72 * t95;
	t88 = t63 * t94;
	t87 = t63 * t93;
	t84 = t69 * t71 * t63 * r_i_i_C(2) + r_i_i_C(1) * t87 + t68 * t89;
	t81 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t73);
	t83 = -pkin(2) + t81 - t96;
	t82 = t81 * qJD(3) - t90;
	t70 = qJ(4) + pkin(7) + pkin(8);
	t78 = -t68 * t56 + r_i_i_C(1) * t88 + qJD(4) * t69 + (t83 * t69 + (-r_i_i_C(3) - t70) * t68) * t72 + t97 * r_i_i_C(2);
	t77 = r_i_i_C(2) * t87 + r_i_i_C(3) * t92 + t68 * qJD(4) + t83 * t93 + (t70 * t72 + t98) * t69;
	t48 = r_i_i_C(2) * t88;
	t1 = [-cos(qJ(1)) * t91 + t78, t78, -t59 * t93 + t69 * t82 + t84, t92, -t69 * t90 + t84; -sin(qJ(1)) * t91 + t77, t77, t48 + t82 * t68 + (t59 + t80) * t92, t93, -t97 * r_i_i_C(1) - t69 * t89 + t48; 0, 0, t98, 0, t79;];
	JaD_transl = t1;
end