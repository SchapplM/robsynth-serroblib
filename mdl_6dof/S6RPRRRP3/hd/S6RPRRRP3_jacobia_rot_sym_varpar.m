% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (218->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(10);
	t35 = cos(t36);
	t54 = t35 ^ 2;
	t43 = cos(qJ(3));
	t34 = sin(t36);
	t41 = sin(qJ(3));
	t49 = t34 * t41;
	t32 = atan2(-t49, -t43);
	t30 = sin(t32);
	t31 = cos(t32);
	t23 = -t30 * t49 - t31 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t53 = t22 * t41;
	t40 = sin(qJ(4));
	t42 = cos(qJ(4));
	t45 = t42 * t43;
	t29 = t34 * t40 + t35 * t45;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t43;
	t28 = -t34 * t42 + t35 * t46;
	t52 = t27 * t28;
	t51 = t30 * t43;
	t37 = t41 ^ 2;
	t47 = t37 / t43 ^ 2;
	t33 = 0.1e1 / (t34 ^ 2 * t47 + 0.1e1);
	t50 = t34 * t33;
	t48 = t35 * t41;
	t44 = t28 ^ 2 * t27 + 0.1e1;
	t38 = 0.1e1 / t43;
	t26 = 0.1e1 / t29;
	t25 = (0.1e1 + t47) * t50;
	t24 = 0.1e1 / t44;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t54 * t37 * t22 + 0.1e1);
	t1 = [t38 * t33 * t48, 0, t25, 0, 0, 0; (-t21 * t49 - (-t31 * t37 * t38 * t50 + (t33 - 0.1e1) * t41 * t30) * t54 * t53) * t20, 0, (t43 * t21 - (-t34 * t51 + t31 * t41 + (-t31 * t49 + t51) * t25) * t53) * t35 * t20, 0, 0, 0; ((-t34 * t46 - t35 * t42) * t26 - (-t34 * t45 + t35 * t40) * t52) * t24, 0, (-t26 * t40 + t42 * t52) * t24 * t48, t44 * t24, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (304->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
	t54 = qJ(1) + pkin(10);
	t51 = cos(t54);
	t71 = t51 ^ 2;
	t60 = cos(qJ(3));
	t50 = sin(t54);
	t59 = sin(qJ(3));
	t66 = t50 * t59;
	t48 = atan2(-t66, -t60);
	t46 = sin(t48);
	t47 = cos(t48);
	t40 = -t46 * t66 - t47 * t60;
	t39 = 0.1e1 / t40 ^ 2;
	t70 = t39 * t59;
	t58 = qJ(4) + qJ(5);
	t52 = sin(t58);
	t53 = cos(t58);
	t63 = t53 * t60;
	t45 = t50 * t52 + t51 * t63;
	t43 = 0.1e1 / t45 ^ 2;
	t64 = t52 * t60;
	t44 = -t50 * t53 + t51 * t64;
	t69 = t43 * t44;
	t68 = t46 * t60;
	t55 = t59 ^ 2;
	t62 = t55 / t60 ^ 2;
	t49 = 0.1e1 / (t50 ^ 2 * t62 + 0.1e1);
	t67 = t50 * t49;
	t65 = t51 * t59;
	t61 = t44 ^ 2 * t43 + 0.1e1;
	t56 = 0.1e1 / t60;
	t42 = 0.1e1 / t45;
	t41 = (0.1e1 + t62) * t67;
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / t61;
	t36 = 0.1e1 / (t71 * t55 * t39 + 0.1e1);
	t35 = t61 * t37;
	t1 = [t56 * t49 * t65, 0, t41, 0, 0, 0; (-t38 * t66 - (-t47 * t55 * t56 * t67 + (t49 - 0.1e1) * t59 * t46) * t71 * t70) * t36, 0, (t60 * t38 - (-t50 * t68 + t47 * t59 + (-t47 * t66 + t68) * t41) * t70) * t51 * t36, 0, 0, 0; ((-t50 * t64 - t51 * t53) * t42 - (-t50 * t63 + t51 * t52) * t69) * t37, 0, (-t42 * t52 + t53 * t69) * t37 * t65, t35, t35, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (905->29), mult. (689->69), div. (139->11), fcn. (1046->9), ass. (0->42)
	t79 = sin(qJ(3));
	t91 = t79 ^ 2;
	t74 = qJ(1) + pkin(10);
	t70 = sin(t74);
	t71 = cos(t74);
	t78 = qJ(4) + qJ(5);
	t73 = cos(t78);
	t72 = sin(t78);
	t80 = cos(qJ(3));
	t84 = t72 * t80;
	t60 = t70 * t84 + t71 * t73;
	t82 = t79 * t72;
	t57 = atan2(-t60, t82);
	t53 = sin(t57);
	t54 = cos(t57);
	t51 = -t53 * t60 + t54 * t82;
	t50 = 0.1e1 / t51 ^ 2;
	t63 = -t70 * t73 + t71 * t84;
	t90 = t50 * t63;
	t89 = t50 * t63 ^ 2;
	t87 = t54 * t60;
	t68 = 0.1e1 / t72;
	t76 = 0.1e1 / t79;
	t86 = t68 * t76;
	t85 = t71 * t79;
	t83 = t73 * t80;
	t64 = t70 * t72 + t71 * t83;
	t59 = 0.1e1 / t64 ^ 2;
	t81 = t59 * t71 ^ 2 * t91;
	t77 = 0.1e1 / t91;
	t69 = 0.1e1 / t72 ^ 2;
	t62 = t70 * t83 - t71 * t72;
	t58 = 0.1e1 / t64;
	t56 = 0.1e1 / (t60 ^ 2 * t69 * t77 + 0.1e1);
	t55 = 0.1e1 / (0.1e1 + t81);
	t52 = t63 * t59 * t55 * t85;
	t49 = 0.1e1 / t51;
	t48 = (t60 * t68 * t77 * t80 + t70) * t56;
	t47 = 0.1e1 / (0.1e1 + t89);
	t46 = (t60 * t69 * t73 - t62 * t68) * t76 * t56;
	t45 = (t64 * t49 - (t54 * t73 * t79 - t53 * t62 + (-t53 * t82 - t87) * t46) * t90) * t47;
	t1 = [-t63 * t56 * t86, 0, t48, t46, t46, 0; (-t60 * t49 - (-t53 + (t86 * t87 + t53) * t56) * t89) * t47, 0, (t48 * t87 * t90 + (-t49 * t85 - (t54 * t80 + (-t48 + t70) * t79 * t53) * t90) * t72) * t47, t45, t45, 0; (-t59 * t62 * t71 + t58 * t70) * t79 * t55, 0, (-t58 * t71 * t80 - t73 * t81) * t55, -t52, -t52, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end