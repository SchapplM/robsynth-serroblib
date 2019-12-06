% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP4
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
%   Wie in S5PRRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (257->15), mult. (243->34), div. (59->8), fcn. (349->9), ass. (0->26)
	t55 = qJ(2) + qJ(3);
	t54 = cos(t55);
	t53 = sin(t55);
	t56 = sin(pkin(8));
	t61 = t56 * t53;
	t49 = atan2(-t61, -t54);
	t47 = sin(t49);
	t64 = t47 * t54;
	t51 = t53 ^ 2;
	t63 = t51 / t54 ^ 2;
	t57 = cos(pkin(8));
	t62 = t54 * t57;
	t58 = sin(qJ(4));
	t59 = cos(qJ(4));
	t46 = t56 * t58 + t59 * t62;
	t44 = 0.1e1 / t46 ^ 2;
	t45 = -t56 * t59 + t58 * t62;
	t60 = t45 ^ 2 * t44 + 0.1e1;
	t48 = cos(t49);
	t43 = 0.1e1 / t60;
	t42 = (0.1e1 + t63) * t56 / (t56 ^ 2 * t63 + 0.1e1);
	t41 = -t47 * t61 - t48 * t54;
	t40 = 0.1e1 / t41 ^ 2;
	t38 = (-t58 / t46 + t59 * t45 * t44) * t57 * t53 * t43;
	t37 = (t54 / t41 - (-t56 * t64 + t48 * t53 + (-t48 * t61 + t64) * t42) * t53 * t40) * t57 / (t57 ^ 2 * t51 * t40 + 0.1e1);
	t1 = [0, t42, t42, 0, 0; 0, t37, t37, 0, 0; 0, t38, t38, t60 * t43, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (374->22), mult. (538->55), div. (111->11), fcn. (804->9), ass. (0->39)
	t73 = qJ(2) + qJ(3);
	t68 = sin(t73);
	t88 = t68 ^ 2;
	t69 = cos(t73);
	t75 = cos(pkin(8));
	t77 = cos(qJ(4));
	t79 = t75 * t77;
	t74 = sin(pkin(8));
	t76 = sin(qJ(4));
	t82 = t74 * t76;
	t60 = t69 * t82 + t79;
	t83 = t68 * t76;
	t58 = atan2(-t60, t83);
	t54 = sin(t58);
	t55 = cos(t58);
	t53 = -t54 * t60 + t55 * t83;
	t52 = 0.1e1 / t53 ^ 2;
	t80 = t75 * t76;
	t81 = t74 * t77;
	t63 = t69 * t80 - t81;
	t87 = t52 * t63;
	t85 = t55 * t60;
	t84 = t68 * t75;
	t64 = t69 * t79 + t82;
	t59 = 0.1e1 / t64 ^ 2;
	t78 = t75 ^ 2 * t88 * t59;
	t72 = 0.1e1 / t76 ^ 2;
	t71 = 0.1e1 / t76;
	t67 = 0.1e1 / t88;
	t62 = t69 * t81 - t80;
	t57 = 0.1e1 / (t60 ^ 2 * t67 * t72 + 0.1e1);
	t56 = 0.1e1 / (0.1e1 + t78);
	t51 = 0.1e1 / t53;
	t50 = (t60 * t67 * t69 * t71 + t74) * t57;
	t49 = 0.1e1 / (t63 ^ 2 * t52 + 0.1e1);
	t48 = (t60 * t72 * t77 - t62 * t71) / t68 * t57;
	t47 = (-t75 * t69 / t64 - t77 * t78) * t56;
	t46 = (t50 * t85 * t87 + (-t51 * t84 - (t55 * t69 + (-t50 + t74) * t54 * t68) * t87) * t76) * t49;
	t1 = [0, t50, t50, t48, 0; 0, t46, t46, (t64 * t51 - (t55 * t68 * t77 - t54 * t62 + (-t54 * t83 - t85) * t48) * t87) * t49, 0; 0, t47, t47, -t63 * t59 * t56 * t84, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end