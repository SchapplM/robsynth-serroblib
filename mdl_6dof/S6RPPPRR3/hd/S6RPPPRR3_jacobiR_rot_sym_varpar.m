% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
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
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [-t5, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (8->4), div. (0->0), fcn. (16->4), ass. (0->7)
	t26 = cos(qJ(1));
	t25 = sin(qJ(1));
	t24 = cos(pkin(9));
	t23 = sin(pkin(9));
	t22 = t25 * t23 + t26 * t24;
	t21 = t26 * t23 - t25 * t24;
	t1 = [t21, 0, 0, 0, 0, 0; t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t22, 0, 0, 0, 0, 0; -t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->4), mult. (16->8), div. (0->0), fcn. (28->6), ass. (0->9)
	t23 = cos(qJ(1));
	t22 = sin(qJ(1));
	t21 = cos(pkin(9));
	t20 = sin(pkin(9));
	t19 = cos(pkin(10));
	t18 = sin(pkin(10));
	t15 = t23 * t20 - t22 * t21;
	t14 = -t22 * t20 - t23 * t21;
	t1 = [t15 * t19, 0, 0, 0, 0, 0; -t14 * t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t15 * t18, 0, 0, 0, 0, 0; t14 * t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->6), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->14)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t32 = sin(pkin(9));
	t33 = cos(pkin(9));
	t24 = -t38 * t32 - t39 * t33;
	t31 = pkin(10) + qJ(5);
	t29 = sin(t31);
	t37 = t24 * t29;
	t30 = cos(t31);
	t36 = t24 * t30;
	t25 = t39 * t32 - t38 * t33;
	t35 = t25 * t29;
	t34 = t25 * t30;
	t1 = [t34, 0, 0, 0, t37, 0; -t36, 0, 0, 0, t35, 0; 0, 0, 0, 0, -t30, 0; -t35, 0, 0, 0, t36, 0; t37, 0, 0, 0, t34, 0; 0, 0, 0, 0, t29, 0; t24, 0, 0, 0, 0, 0; t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->16), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->20)
	t102 = cos(qJ(1));
	t101 = sin(qJ(1));
	t90 = pkin(10) + qJ(5);
	t88 = sin(t90);
	t91 = sin(qJ(6));
	t100 = t88 * t91;
	t92 = cos(qJ(6));
	t99 = t88 * t92;
	t89 = cos(t90);
	t98 = t89 * t91;
	t97 = t89 * t92;
	t96 = cos(pkin(9));
	t95 = sin(pkin(9));
	t83 = -t101 * t95 - t102 * t96;
	t84 = -t101 * t96 + t102 * t95;
	t94 = t83 * t91 + t84 * t97;
	t93 = -t83 * t92 + t84 * t98;
	t82 = -t83 * t97 + t84 * t91;
	t81 = t83 * t98 + t84 * t92;
	t1 = [t94, 0, 0, 0, t83 * t99, t81; t82, 0, 0, 0, t84 * t99, t93; 0, 0, 0, 0, -t97, t100; -t93, 0, 0, 0, -t83 * t100, -t82; t81, 0, 0, 0, -t84 * t100, t94; 0, 0, 0, 0, t98, t99; t84 * t88, 0, 0, 0, -t83 * t89, 0; -t83 * t88, 0, 0, 0, -t84 * t89, 0; 0, 0, 0, 0, -t88, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end