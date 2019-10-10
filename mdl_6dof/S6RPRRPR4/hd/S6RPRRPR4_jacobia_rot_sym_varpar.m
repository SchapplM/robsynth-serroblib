% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR4
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
%   Wie in S6RPRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
	t63 = pkin(10) + qJ(3) + qJ(4);
	t62 = cos(t63);
	t61 = sin(t63);
	t66 = sin(qJ(1));
	t72 = t66 * t61;
	t56 = atan2(-t72, -t62);
	t50 = sin(t56);
	t51 = cos(t56);
	t47 = -t50 * t72 - t51 * t62;
	t46 = 0.1e1 / t47 ^ 2;
	t67 = cos(qJ(1));
	t78 = t46 * t67 ^ 2;
	t77 = t50 * t62;
	t65 = cos(pkin(11));
	t68 = t67 * t65;
	t64 = sin(pkin(11));
	t71 = t66 * t64;
	t55 = t62 * t68 + t71;
	t53 = 0.1e1 / t55 ^ 2;
	t69 = t67 * t64;
	t70 = t66 * t65;
	t54 = t62 * t69 - t70;
	t76 = t53 * t54;
	t58 = t61 ^ 2;
	t75 = t58 / t62 ^ 2;
	t74 = t61 * t67;
	t57 = 0.1e1 / (t66 ^ 2 * t75 + 0.1e1);
	t73 = t66 * t57;
	t59 = 0.1e1 / t62;
	t52 = 0.1e1 / t55;
	t49 = 0.1e1 / (t54 ^ 2 * t53 + 0.1e1);
	t48 = (0.1e1 + t75) * t73;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t58 * t78 + 0.1e1);
	t43 = (-t52 * t64 + t65 * t76) * t49 * t74;
	t42 = (t62 * t45 - (-t66 * t77 + t51 * t61 + (-t51 * t72 + t77) * t48) * t61 * t46) * t67 * t44;
	t1 = [t59 * t57 * t74, 0, t48, t48, 0, 0; (-t45 * t72 - (-t51 * t58 * t59 * t73 + (t57 - 0.1e1) * t61 * t50) * t61 * t78) * t44, 0, t42, t42, 0, 0; ((-t62 * t71 - t68) * t52 - (-t62 * t70 + t69) * t76) * t49, 0, t43, t43, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (618->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t70 = pkin(10) + qJ(3) + qJ(4);
	t67 = cos(t70);
	t66 = sin(t70);
	t72 = sin(qJ(1));
	t79 = t72 * t66;
	t61 = atan2(-t79, -t67);
	t59 = sin(t61);
	t60 = cos(t61);
	t52 = -t59 * t79 - t60 * t67;
	t51 = 0.1e1 / t52 ^ 2;
	t73 = cos(qJ(1));
	t85 = t51 * t73 ^ 2;
	t71 = pkin(11) + qJ(6);
	t69 = cos(t71);
	t75 = t73 * t69;
	t68 = sin(t71);
	t78 = t72 * t68;
	t58 = t67 * t75 + t78;
	t56 = 0.1e1 / t58 ^ 2;
	t76 = t73 * t68;
	t77 = t72 * t69;
	t57 = t67 * t76 - t77;
	t84 = t56 * t57;
	t83 = t59 * t67;
	t63 = t66 ^ 2;
	t82 = t63 / t67 ^ 2;
	t81 = t66 * t73;
	t62 = 0.1e1 / (t72 ^ 2 * t82 + 0.1e1);
	t80 = t72 * t62;
	t74 = t57 ^ 2 * t56 + 0.1e1;
	t64 = 0.1e1 / t67;
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t74;
	t53 = (0.1e1 + t82) * t80;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (t63 * t85 + 0.1e1);
	t48 = (-t55 * t68 + t69 * t84) * t54 * t81;
	t47 = (t67 * t50 - (-t72 * t83 + t60 * t66 + (-t60 * t79 + t83) * t53) * t66 * t51) * t73 * t49;
	t1 = [t64 * t62 * t81, 0, t53, t53, 0, 0; (-t50 * t79 - (-t60 * t63 * t64 * t80 + (t62 - 0.1e1) * t66 * t59) * t66 * t85) * t49, 0, t47, t47, 0, 0; ((-t67 * t78 - t75) * t55 - (-t67 * t77 + t76) * t84) * t54, 0, t48, t48, 0, t74 * t54;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end