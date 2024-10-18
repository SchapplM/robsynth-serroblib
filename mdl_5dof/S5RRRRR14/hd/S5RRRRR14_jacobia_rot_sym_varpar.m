% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR14
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
%   Wie in S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRR14_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (164->15), mult. (205->35), div. (45->9), fcn. (315->9), ass. (0->27)
	t47 = cos(pkin(5));
	t45 = qJ(1) + qJ(2);
	t42 = cos(t45);
	t46 = sin(pkin(5));
	t54 = t42 * t46;
	t40 = atan2(t54, t47);
	t37 = sin(t40);
	t38 = cos(t40);
	t33 = t37 * t54 + t38 * t47;
	t41 = sin(t45);
	t56 = 0.1e1 / t33 ^ 2 * t41 ^ 2;
	t43 = t46 ^ 2;
	t39 = 0.1e1 / (0.1e1 + t42 ^ 2 * t43 / t47 ^ 2);
	t55 = t39 / t47;
	t48 = sin(qJ(3));
	t53 = t47 * t48;
	t49 = cos(qJ(3));
	t52 = t47 * t49;
	t36 = -t41 * t53 + t42 * t49;
	t34 = 0.1e1 / t36 ^ 2;
	t35 = t41 * t52 + t42 * t48;
	t51 = t35 ^ 2 * t34 + 0.1e1;
	t50 = t41 * t46 * t55;
	t32 = 0.1e1 / t51;
	t29 = ((-t41 * t48 + t42 * t52) / t36 - (-t41 * t49 - t42 * t53) * t35 * t34) * t32;
	t28 = (0.1e1 / t33 * t54 - (-t38 * t42 * t43 * t55 + (t39 - 0.1e1) * t46 * t37) * t46 * t56) / (t43 * t56 + 0.1e1);
	t1 = [-t50, -t50, 0, 0, 0; t28, t28, 0, 0, 0; t29, t29, t51 * t32, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (416->21), mult. (280->38), div. (52->9), fcn. (370->13), ass. (0->32)
	t78 = cos(pkin(5));
	t76 = qJ(1) + qJ(2);
	t72 = cos(t76);
	t77 = sin(pkin(5));
	t80 = t72 * t77;
	t65 = atan2(t80, t78);
	t60 = sin(t65);
	t61 = cos(t65);
	t55 = t60 * t80 + t61 * t78;
	t70 = sin(t76);
	t83 = 0.1e1 / t55 ^ 2 * t70 ^ 2;
	t75 = qJ(3) + qJ(4);
	t67 = pkin(5) + t75;
	t68 = pkin(5) - t75;
	t62 = sin(t67) / 0.2e1 - sin(t68) / 0.2e1;
	t71 = cos(t75);
	t59 = -t70 * t62 + t72 * t71;
	t57 = 0.1e1 / t59 ^ 2;
	t63 = cos(t68) / 0.2e1 + cos(t67) / 0.2e1;
	t69 = sin(t75);
	t58 = t70 * t63 + t72 * t69;
	t82 = t58 ^ 2 * t57;
	t73 = t77 ^ 2;
	t64 = 0.1e1 / (0.1e1 + t72 ^ 2 * t73 / t78 ^ 2);
	t81 = t64 / t78;
	t79 = t70 * t77 * t81;
	t56 = 0.1e1 / t59;
	t52 = 0.1e1 / (0.1e1 + t82);
	t51 = (t59 * t56 + t82) * t52;
	t50 = ((t72 * t63 - t70 * t69) * t56 - (-t72 * t62 - t70 * t71) * t58 * t57) * t52;
	t49 = (0.1e1 / t55 * t80 - (-t61 * t72 * t73 * t81 + (t64 - 0.1e1) * t77 * t60) * t77 * t83) / (t73 * t83 + 0.1e1);
	t1 = [-t79, -t79, 0, 0, 0; t49, t49, 0, 0, 0; t50, t50, t51, t51, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (635->21), mult. (322->38), div. (58->9), fcn. (415->13), ass. (0->32)
	t88 = cos(pkin(5));
	t86 = qJ(1) + qJ(2);
	t82 = cos(t86);
	t87 = sin(pkin(5));
	t90 = t82 * t87;
	t75 = atan2(t90, t88);
	t72 = sin(t75);
	t73 = cos(t75);
	t65 = t72 * t90 + t73 * t88;
	t81 = sin(t86);
	t93 = 0.1e1 / t65 ^ 2 * t81 ^ 2;
	t83 = qJ(3) + qJ(4) + qJ(5);
	t77 = pkin(5) + t83;
	t78 = pkin(5) - t83;
	t70 = sin(t77) / 0.2e1 - sin(t78) / 0.2e1;
	t80 = cos(t83);
	t69 = -t81 * t70 + t82 * t80;
	t67 = 0.1e1 / t69 ^ 2;
	t71 = cos(t78) / 0.2e1 + cos(t77) / 0.2e1;
	t79 = sin(t83);
	t68 = t81 * t71 + t82 * t79;
	t92 = t68 ^ 2 * t67;
	t84 = t87 ^ 2;
	t74 = 0.1e1 / (0.1e1 + t82 ^ 2 * t84 / t88 ^ 2);
	t91 = t74 / t88;
	t89 = t81 * t87 * t91;
	t66 = 0.1e1 / t69;
	t62 = 0.1e1 / (0.1e1 + t92);
	t61 = (0.1e1 / t65 * t90 - (-t73 * t82 * t84 * t91 + (t74 - 0.1e1) * t87 * t72) * t87 * t93) / (t84 * t93 + 0.1e1);
	t60 = (t69 * t66 + t92) * t62;
	t59 = ((t82 * t71 - t81 * t79) * t66 - (-t82 * t70 - t81 * t80) * t68 * t67) * t62;
	t1 = [-t89, -t89, 0, 0, 0; t61, t61, 0, 0, 0; t59, t59, t60, t60, t60;];
	Ja_rot = t1;
end