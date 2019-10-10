% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP2
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
%   Wie in S6RPRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t49 = qJ(3) + pkin(10);
	t47 = cos(t49);
	t45 = sin(t49);
	t50 = qJ(1) + pkin(9);
	t46 = sin(t50);
	t58 = t46 * t45;
	t40 = atan2(-t58, -t47);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t58 - t39 * t47;
	t30 = 0.1e1 / t31 ^ 2;
	t48 = cos(t50);
	t64 = t30 * t48 ^ 2;
	t52 = cos(qJ(5));
	t54 = t48 * t52;
	t51 = sin(qJ(5));
	t57 = t46 * t51;
	t37 = t47 * t54 + t57;
	t35 = 0.1e1 / t37 ^ 2;
	t55 = t48 * t51;
	t56 = t46 * t52;
	t36 = t47 * t55 - t56;
	t63 = t35 * t36;
	t62 = t38 * t47;
	t42 = t45 ^ 2;
	t61 = t42 / t47 ^ 2;
	t60 = t45 * t48;
	t41 = 0.1e1 / (t46 ^ 2 * t61 + 0.1e1);
	t59 = t46 * t41;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t43 = 0.1e1 / t47;
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t53;
	t32 = (0.1e1 + t61) * t59;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t42 * t64 + 0.1e1);
	t1 = [t43 * t41 * t60, 0, t32, 0, 0, 0; (-t29 * t58 - (-t39 * t42 * t43 * t59 + (t41 - 0.1e1) * t45 * t38) * t45 * t64) * t28, 0, (t47 * t29 - (-t46 * t62 + t39 * t45 + (-t39 * t58 + t62) * t32) * t45 * t30) * t48 * t28, 0, 0, 0; ((-t47 * t57 - t54) * t34 - (-t47 * t56 + t55) * t63) * t33, 0, (-t34 * t51 + t52 * t63) * t33 * t60, 0, t53 * t33, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (572->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
	t60 = qJ(3) + pkin(10);
	t56 = sin(t60);
	t78 = t56 ^ 2;
	t58 = cos(t60);
	t61 = qJ(1) + pkin(9);
	t59 = cos(t61);
	t65 = cos(qJ(5));
	t67 = t59 * t65;
	t57 = sin(t61);
	t64 = sin(qJ(5));
	t70 = t57 * t64;
	t45 = t58 * t70 + t67;
	t71 = t56 * t64;
	t42 = atan2(-t45, t71);
	t39 = sin(t42);
	t40 = cos(t42);
	t37 = -t39 * t45 + t40 * t71;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t59 * t64;
	t69 = t57 * t65;
	t48 = t58 * t68 - t69;
	t77 = t36 * t48;
	t75 = t40 * t45;
	t74 = t48 ^ 2 * t36;
	t53 = 0.1e1 / t56;
	t62 = 0.1e1 / t64;
	t73 = t53 * t62;
	t72 = t56 * t59;
	t49 = t58 * t67 + t70;
	t44 = 0.1e1 / t49 ^ 2;
	t66 = t59 ^ 2 * t78 * t44;
	t63 = 0.1e1 / t64 ^ 2;
	t54 = 0.1e1 / t78;
	t47 = t58 * t69 - t68;
	t43 = 0.1e1 / t49;
	t41 = 0.1e1 / (t45 ^ 2 * t54 * t63 + 0.1e1);
	t38 = 0.1e1 / (0.1e1 + t66);
	t35 = 0.1e1 / t37;
	t34 = (t45 * t54 * t58 * t62 + t57) * t41;
	t33 = 0.1e1 / (0.1e1 + t74);
	t32 = (t45 * t63 * t65 - t47 * t62) * t53 * t41;
	t1 = [-t48 * t41 * t73, 0, t34, 0, t32, 0; (-t45 * t35 - (-t39 + (t73 * t75 + t39) * t41) * t74) * t33, 0, (t34 * t75 * t77 + (-t35 * t72 - (t40 * t58 + (-t34 + t57) * t56 * t39) * t77) * t64) * t33, 0, (t49 * t35 - (t40 * t56 * t65 - t39 * t47 + (-t39 * t71 - t75) * t32) * t77) * t33, 0; (-t44 * t47 * t59 + t43 * t57) * t56 * t38, 0, (-t43 * t58 * t59 - t65 * t66) * t38, 0, -t48 * t44 * t38 * t72, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end