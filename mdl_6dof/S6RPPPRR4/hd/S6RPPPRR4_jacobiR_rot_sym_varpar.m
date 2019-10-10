% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (8->4), div. (0->0), fcn. (16->4), ass. (0->7)
	t18 = cos(qJ(1));
	t17 = sin(qJ(1));
	t16 = cos(pkin(9));
	t15 = sin(pkin(9));
	t14 = t18 * t15 - t17 * t16;
	t13 = -t17 * t15 - t18 * t16;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->5), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t28 = sin(pkin(9));
	t29 = cos(pkin(9));
	t22 = -t33 * t28 - t34 * t29;
	t26 = sin(qJ(5));
	t32 = t22 * t26;
	t27 = cos(qJ(5));
	t31 = t22 * t27;
	t23 = t34 * t28 - t33 * t29;
	t30 = t23 * t26;
	t21 = t23 * t27;
	t1 = [t32, 0, 0, 0, t21, 0; t30, 0, 0, 0, -t31, 0; 0, 0, 0, 0, t26, 0; t31, 0, 0, 0, -t30, 0; t21, 0, 0, 0, t32, 0; 0, 0, 0, 0, t27, 0; t23, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (38->15), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
	t94 = cos(qJ(1));
	t93 = sin(qJ(1));
	t81 = sin(qJ(6));
	t82 = sin(qJ(5));
	t92 = t82 * t81;
	t83 = cos(qJ(6));
	t91 = t82 * t83;
	t84 = cos(qJ(5));
	t90 = t84 * t81;
	t89 = t84 * t83;
	t88 = cos(pkin(9));
	t87 = sin(pkin(9));
	t76 = -t87 * t93 - t88 * t94;
	t77 = t87 * t94 - t88 * t93;
	t86 = t76 * t91 + t77 * t81;
	t85 = t76 * t92 - t77 * t83;
	t75 = -t76 * t81 + t77 * t91;
	t74 = -t76 * t83 - t77 * t92;
	t1 = [t86, 0, 0, 0, t77 * t89, t74; t75, 0, 0, 0, -t76 * t89, t85; 0, 0, 0, 0, t91, t90; -t85, 0, 0, 0, -t77 * t90, -t75; t74, 0, 0, 0, t76 * t90, t86; 0, 0, 0, 0, -t92, t89; -t76 * t84, 0, 0, 0, t77 * t82, 0; -t77 * t84, 0, 0, 0, -t76 * t82, 0; 0, 0, 0, 0, -t84, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end