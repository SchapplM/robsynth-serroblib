% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(9);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t23 = t18 * t15;
	t16 = cos(t17);
	t22 = t18 * t16;
	t19 = cos(qJ(1));
	t21 = t19 * t15;
	t20 = t19 * t16;
	t1 = [-t22, -t21, 0, 0, 0, 0; t20, -t23, 0, 0, 0, 0; 0, t16, 0, 0, 0, 0; t23, -t20, 0, 0, 0, 0; -t21, -t22, 0, 0, 0, 0; 0, -t15, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t50 = qJ(2) + pkin(9);
	t48 = sin(t50);
	t51 = sin(qJ(1));
	t54 = t51 * t48;
	t49 = cos(t50);
	t52 = cos(qJ(1));
	t53 = t52 * t49;
	t47 = t52 * t48;
	t46 = t51 * t49;
	t1 = [t52, 0, 0, 0, 0, 0; t51, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t46, t47, 0, 0, 0, 0; -t53, t54, 0, 0, 0, 0; 0, -t49, 0, 0, 0, 0; -t54, t53, 0, 0, 0, 0; t47, t46, 0, 0, 0, 0; 0, t48, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t57 = sin(pkin(10));
	t59 = sin(qJ(1));
	t64 = t59 * t57;
	t58 = cos(pkin(10));
	t63 = t59 * t58;
	t60 = cos(qJ(1));
	t62 = t60 * t57;
	t61 = t60 * t58;
	t56 = qJ(2) + pkin(9);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [-t54 * t64 + t61, t55 * t62, 0, 0, 0, 0; t54 * t62 + t63, t55 * t64, 0, 0, 0, 0; 0, t54 * t57, 0, 0, 0, 0; -t54 * t63 - t62, t55 * t61, 0, 0, 0, 0; t54 * t61 - t64, t55 * t63, 0, 0, 0, 0; 0, t54 * t58, 0, 0, 0, 0; -t59 * t55, -t60 * t54, 0, 0, 0, 0; t60 * t55, -t59 * t54, 0, 0, 0, 0; 0, t55, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (57->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t83 = qJ(2) + pkin(9);
	t79 = sin(t83);
	t84 = sin(qJ(1));
	t91 = t84 * t79;
	t82 = pkin(10) + qJ(6);
	t80 = cos(t82);
	t90 = t84 * t80;
	t81 = cos(t83);
	t89 = t84 * t81;
	t85 = cos(qJ(1));
	t88 = t85 * t79;
	t87 = t85 * t80;
	t86 = t85 * t81;
	t78 = sin(t82);
	t77 = -t78 * t91 + t87;
	t76 = t85 * t78 + t79 * t90;
	t75 = t78 * t88 + t90;
	t74 = -t84 * t78 + t79 * t87;
	t1 = [t77, t78 * t86, 0, 0, 0, t74; t75, t78 * t89, 0, 0, 0, t76; 0, t79 * t78, 0, 0, 0, -t81 * t80; -t76, t80 * t86, 0, 0, 0, -t75; t74, t80 * t89, 0, 0, 0, t77; 0, t79 * t80, 0, 0, 0, t81 * t78; -t89, -t88, 0, 0, 0, 0; t86, -t91, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end