% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:26
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:45
	% EndTime: 2019-10-24 10:26:45
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (222->18), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
	t56 = sin(pkin(9));
	t57 = sin(pkin(5));
	t68 = t56 * t57;
	t60 = cos(pkin(5));
	t55 = sin(pkin(10));
	t58 = cos(pkin(10));
	t62 = sin(qJ(2));
	t64 = cos(qJ(2));
	t66 = t64 * t55 + t62 * t58;
	t52 = t66 * t60;
	t53 = t62 * t55 - t64 * t58;
	t59 = cos(pkin(9));
	t47 = -t56 * t52 - t59 * t53;
	t61 = sin(qJ(4));
	t63 = cos(qJ(4));
	t42 = t47 * t63 + t61 * t68;
	t40 = 0.1e1 / t42 ^ 2;
	t41 = t47 * t61 - t63 * t68;
	t67 = t41 ^ 2 * t40 + 0.1e1;
	t65 = t53 * t60;
	t51 = t66 * t57;
	t50 = t53 * t57;
	t49 = 0.1e1 / t50 ^ 2;
	t45 = t56 * t65 - t59 * t66;
	t44 = -t56 * t66 - t59 * t65;
	t43 = -t59 * t52 + t56 * t53;
	t39 = atan2(t44, t50);
	t37 = cos(t39);
	t36 = sin(t39);
	t35 = 0.1e1 / t67;
	t34 = t36 * t44 + t37 * t50;
	t33 = 0.1e1 / t34 ^ 2;
	t31 = (t43 / t50 - t51 * t44 * t49) / (t44 ^ 2 * t49 + 0.1e1);
	t1 = [0, t31, 0, 0, 0; 0, (t47 / t34 + (t36 * t43 + t37 * t51 + (-t36 * t50 + t37 * t44) * t31) * t45 * t33) / (t45 ^ 2 * t33 + 0.1e1), 0, 0, 0; 0, (t61 / t42 - t63 * t41 * t40) * t45 * t35, 0, t67 * t35, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:45
	% EndTime: 2019-10-24 10:26:45
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (632->32), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->54)
	t85 = sin(pkin(5));
	t93 = cos(qJ(4));
	t100 = t85 * t93;
	t88 = cos(pkin(5));
	t83 = sin(pkin(10));
	t86 = cos(pkin(10));
	t91 = sin(qJ(2));
	t94 = cos(qJ(2));
	t96 = t94 * t83 + t91 * t86;
	t78 = t96 * t88;
	t79 = t91 * t83 - t94 * t86;
	t84 = sin(pkin(9));
	t87 = cos(pkin(9));
	t67 = t87 * t78 - t84 * t79;
	t90 = sin(qJ(4));
	t61 = t87 * t100 + t67 * t90;
	t77 = t96 * t85;
	t73 = t77 * t90 - t88 * t93;
	t60 = atan2(-t61, t73);
	t57 = sin(t60);
	t58 = cos(t60);
	t51 = -t57 * t61 + t58 * t73;
	t50 = 0.1e1 / t51 ^ 2;
	t97 = -t84 * t78 - t87 * t79;
	t64 = -t84 * t100 + t90 * t97;
	t105 = t50 * t64;
	t101 = t85 * t90;
	t65 = t84 * t101 + t93 * t97;
	t95 = t79 * t88;
	t69 = t84 * t95 - t87 * t96;
	t89 = sin(qJ(5));
	t92 = cos(qJ(5));
	t56 = t65 * t92 - t69 * t89;
	t54 = 0.1e1 / t56 ^ 2;
	t55 = t65 * t89 + t69 * t92;
	t104 = t54 * t55;
	t72 = 0.1e1 / t73 ^ 2;
	t103 = t61 * t72;
	t102 = t69 * t93;
	t99 = t55 ^ 2 * t54 + 0.1e1;
	t98 = -t57 * t73 - t58 * t61;
	t76 = t79 * t85;
	t74 = t77 * t93 + t88 * t90;
	t71 = 0.1e1 / t73;
	t66 = -t84 * t96 - t87 * t95;
	t63 = -t87 * t101 + t67 * t93;
	t59 = 0.1e1 / (t61 ^ 2 * t72 + 0.1e1);
	t53 = 0.1e1 / t56;
	t52 = 0.1e1 / t99;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (t64 ^ 2 * t50 + 0.1e1);
	t47 = (-t76 * t103 - t66 * t71) * t90 * t59;
	t46 = (t74 * t103 - t63 * t71) * t59;
	t1 = [0, t47, 0, t46, 0; 0, (t69 * t90 * t49 - ((-t57 * t66 - t58 * t76) * t90 + t98 * t47) * t105) * t48, 0, (t65 * t49 - (t98 * t46 - t57 * t63 + t58 * t74) * t105) * t48, 0; 0, ((t89 * t102 - t92 * t97) * t53 - (t92 * t102 + t89 * t97) * t104) * t52, 0, (t92 * t104 - t53 * t89) * t64 * t52, t99 * t52;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end