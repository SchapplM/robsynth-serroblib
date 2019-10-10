% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP8
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
%   Wie in S6RPRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (323->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t63 = qJ(3) + qJ(4);
	t61 = sin(t63);
	t62 = cos(t63);
	t67 = cos(qJ(1));
	t71 = t67 * t62;
	t56 = atan2(-t71, t61);
	t54 = sin(t56);
	t55 = cos(t56);
	t47 = -t54 * t71 + t55 * t61;
	t46 = 0.1e1 / t47 ^ 2;
	t65 = sin(qJ(1));
	t79 = t46 * t65 ^ 2;
	t64 = sin(qJ(5));
	t70 = t67 * t64;
	t66 = cos(qJ(5));
	t73 = t65 * t66;
	t53 = t61 * t73 + t70;
	t51 = 0.1e1 / t53 ^ 2;
	t69 = t67 * t66;
	t74 = t65 * t64;
	t52 = t61 * t74 - t69;
	t78 = t51 * t52;
	t77 = t54 * t61;
	t60 = t62 ^ 2;
	t76 = 0.1e1 / t61 ^ 2 * t60;
	t75 = t62 * t65;
	t57 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
	t72 = t67 * t57;
	t68 = t52 ^ 2 * t51 + 0.1e1;
	t58 = 0.1e1 / t61;
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t68;
	t48 = (0.1e1 + t76) * t72;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t60 * t79 + 0.1e1);
	t43 = (t50 * t64 - t66 * t78) * t49 * t75;
	t42 = (t61 * t45 + (t67 * t77 + t55 * t62 + (-t55 * t71 - t77) * t48) * t62 * t46) * t65 * t44;
	t1 = [t58 * t57 * t75, 0, t48, t48, 0, 0; (-t45 * t71 + (-t55 * t58 * t60 * t72 + (-t57 + 0.1e1) * t62 * t54) * t62 * t79) * t44, 0, t42, t42, 0, 0; ((t61 * t70 + t73) * t50 - (t61 * t69 - t74) * t78) * t49, 0, t43, t43, t68 * t49, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (478->25), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
	t74 = qJ(3) + qJ(4);
	t70 = cos(t74);
	t91 = t70 ^ 2;
	t69 = sin(t74);
	t75 = sin(qJ(5));
	t78 = cos(qJ(1));
	t81 = t78 * t75;
	t76 = sin(qJ(1));
	t77 = cos(qJ(5));
	t82 = t76 * t77;
	t63 = t69 * t81 + t82;
	t85 = t70 * t75;
	t58 = atan2(t63, t85);
	t54 = sin(t58);
	t55 = cos(t58);
	t53 = t54 * t63 + t55 * t85;
	t52 = 0.1e1 / t53 ^ 2;
	t80 = t78 * t77;
	t83 = t76 * t75;
	t61 = t69 * t83 - t80;
	t90 = t52 * t61;
	t88 = t55 * t63;
	t87 = t61 ^ 2 * t52;
	t67 = 0.1e1 / t70;
	t71 = 0.1e1 / t75;
	t86 = t67 * t71;
	t84 = t70 * t76;
	t62 = t69 * t82 + t81;
	t60 = 0.1e1 / t62 ^ 2;
	t79 = t76 ^ 2 * t91 * t60;
	t72 = 0.1e1 / t75 ^ 2;
	t68 = 0.1e1 / t91;
	t64 = t69 * t80 - t83;
	t59 = 0.1e1 / t62;
	t57 = 0.1e1 / (t63 ^ 2 * t68 * t72 + 0.1e1);
	t56 = 0.1e1 / (0.1e1 + t79);
	t51 = 0.1e1 / t53;
	t50 = (t63 * t68 * t69 * t71 + t78) * t57;
	t49 = 0.1e1 / (0.1e1 + t87);
	t48 = (-t63 * t72 * t77 + t64 * t71) * t67 * t57;
	t47 = (-t59 * t69 * t76 - t77 * t79) * t56;
	t46 = (-t50 * t88 * t90 + (t51 * t84 - (-t55 * t69 + (-t50 + t78) * t54 * t70) * t90) * t75) * t49;
	t1 = [-t61 * t57 * t86, 0, t50, t50, t48, 0; (t63 * t51 - (-t54 + (-t86 * t88 + t54) * t57) * t87) * t49, 0, t46, t46, (t62 * t51 - (t55 * t70 * t77 + t54 * t64 + (-t54 * t85 + t88) * t48) * t90) * t49, 0; (-t60 * t64 * t76 + t59 * t78) * t70 * t56, 0, t47, t47, t61 * t60 * t56 * t84, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end