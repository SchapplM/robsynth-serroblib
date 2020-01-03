% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t11 = qJ(1) + qJ(2);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [0, 0, 0, 0, 0; t10, t10, 0, 0, 0; t9, t9, 0, 0, 0; 0, 0, 0, 0, 0; -t9, -t9, 0, 0, 0; t10, t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t29 = qJ(1) + qJ(2);
	t27 = sin(t29);
	t30 = sin(pkin(9));
	t33 = t27 * t30;
	t28 = cos(t29);
	t32 = t28 * t30;
	t31 = cos(pkin(9));
	t26 = t28 * t31;
	t25 = t27 * t31;
	t1 = [0, 0, 0, 0, 0; t26, t26, 0, 0, 0; t25, t25, 0, 0, 0; 0, 0, 0, 0, 0; -t32, -t32, 0, 0, 0; -t33, -t33, 0, 0, 0; 0, 0, 0, 0, 0; t27, t27, 0, 0, 0; -t28, -t28, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->10), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t68 = cos(pkin(9));
	t69 = sin(qJ(4));
	t72 = t68 * t69;
	t70 = cos(qJ(4));
	t71 = t68 * t70;
	t67 = sin(pkin(9));
	t66 = qJ(1) + qJ(2);
	t65 = cos(t66);
	t64 = sin(t66);
	t62 = t65 * t67;
	t61 = t64 * t67;
	t60 = t64 * t69 + t65 * t71;
	t59 = -t64 * t70 + t65 * t72;
	t58 = t64 * t71 - t65 * t69;
	t57 = -t64 * t72 - t65 * t70;
	t1 = [0, 0, 0, -t67 * t69, 0; t60, t60, 0, t57, 0; t58, t58, 0, t59, 0; 0, 0, 0, -t67 * t70, 0; -t59, -t59, 0, -t58, 0; t57, t57, 0, t60, 0; 0, 0, 0, 0, 0; t62, t62, 0, 0, 0; t61, t61, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (92->14), mult. (56->14), div. (0->0), fcn. (96->6), ass. (0->19)
	t87 = qJ(1) + qJ(2);
	t83 = sin(t87);
	t89 = cos(pkin(9));
	t93 = t83 * t89;
	t85 = cos(t87);
	t92 = t85 * t89;
	t86 = qJ(4) + qJ(5);
	t82 = sin(t86);
	t88 = sin(pkin(9));
	t91 = t88 * t82;
	t84 = cos(t86);
	t90 = t88 * t84;
	t81 = t85 * t88;
	t80 = t83 * t88;
	t77 = t83 * t82 + t84 * t92;
	t76 = t82 * t92 - t83 * t84;
	t75 = -t85 * t82 + t84 * t93;
	t74 = -t82 * t93 - t85 * t84;
	t1 = [0, 0, 0, -t91, -t91; t77, t77, 0, t74, t74; t75, t75, 0, t76, t76; 0, 0, 0, -t90, -t90; -t76, -t76, 0, -t75, -t75; t74, t74, 0, t77, t77; 0, 0, 0, 0, 0; t81, t81, 0, 0, 0; t80, t80, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end