% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
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
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
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
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t38 = sin(qJ(2));
	t39 = sin(qJ(1));
	t44 = t39 * t38;
	t40 = cos(qJ(2));
	t43 = t39 * t40;
	t41 = cos(qJ(1));
	t42 = t41 * t38;
	t37 = t41 * t40;
	t1 = [-t43, -t42, 0, 0, 0, 0; t37, -t44, 0, 0, 0, 0; 0, t40, 0, 0, 0, 0; t41, 0, 0, 0, 0, 0; t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t44, t37, 0, 0, 0, 0; t42, t43, 0, 0, 0, 0; 0, t38, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->12), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->13)
	t35 = sin(qJ(4));
	t36 = sin(qJ(2));
	t38 = cos(qJ(4));
	t39 = cos(qJ(2));
	t42 = t39 * t35 - t36 * t38;
	t41 = t36 * t35 + t39 * t38;
	t40 = cos(qJ(1));
	t37 = sin(qJ(1));
	t31 = t41 * t40;
	t30 = t42 * t40;
	t29 = t41 * t37;
	t28 = t42 * t37;
	t1 = [-t29, t30, 0, -t30, 0, 0; t31, t28, 0, -t28, 0, 0; 0, t41, 0, -t41, 0, 0; t28, t31, 0, -t31, 0, 0; -t30, t29, 0, -t29, 0, 0; 0, -t42, 0, t42, 0, 0; -t40, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->20), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->25)
	t115 = sin(qJ(4));
	t118 = cos(qJ(2));
	t132 = sin(qJ(2));
	t133 = cos(qJ(4));
	t109 = -t118 * t115 + t132 * t133;
	t108 = t132 * t115 + t118 * t133;
	t116 = sin(qJ(1));
	t104 = t109 * t116;
	t114 = sin(qJ(5));
	t131 = t104 * t114;
	t117 = cos(qJ(5));
	t130 = t104 * t117;
	t119 = cos(qJ(1));
	t107 = t109 * t119;
	t129 = t107 * t114;
	t128 = t107 * t117;
	t127 = t108 * t114;
	t126 = t108 * t117;
	t105 = t108 * t116;
	t121 = -t105 * t117 - t119 * t114;
	t120 = t105 * t114 - t119 * t117;
	t106 = t108 * t119;
	t103 = t106 * t117 - t116 * t114;
	t102 = -t106 * t114 - t116 * t117;
	t1 = [t121, -t128, 0, t128, t102, 0; t103, -t130, 0, t130, -t120, 0; 0, t126, 0, -t126, -t109 * t114, 0; t120, t129, 0, -t129, -t103, 0; t102, t131, 0, -t131, t121, 0; 0, -t127, 0, t127, -t109 * t117, 0; t104, -t106, 0, t106, 0, 0; -t107, -t105, 0, t105, 0, 0; 0, -t109, 0, t109, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (98->25), mult. (164->24), div. (0->0), fcn. (248->8), ass. (0->28)
	t138 = sin(qJ(4));
	t140 = cos(qJ(2));
	t155 = sin(qJ(2));
	t156 = cos(qJ(4));
	t129 = -t140 * t138 + t155 * t156;
	t128 = t155 * t138 + t140 * t156;
	t139 = sin(qJ(1));
	t124 = t129 * t139;
	t137 = qJ(5) + qJ(6);
	t135 = sin(t137);
	t154 = t124 * t135;
	t136 = cos(t137);
	t153 = t124 * t136;
	t141 = cos(qJ(1));
	t127 = t129 * t141;
	t152 = t127 * t135;
	t151 = t127 * t136;
	t150 = t128 * t135;
	t149 = t128 * t136;
	t148 = t129 * t135;
	t147 = t129 * t136;
	t125 = t128 * t139;
	t121 = -t125 * t136 - t141 * t135;
	t142 = t125 * t135 - t141 * t136;
	t126 = t128 * t141;
	t123 = t126 * t136 - t139 * t135;
	t122 = -t126 * t135 - t139 * t136;
	t1 = [t121, -t151, 0, t151, t122, t122; t123, -t153, 0, t153, -t142, -t142; 0, t149, 0, -t149, -t148, -t148; t142, t152, 0, -t152, -t123, -t123; t122, t154, 0, -t154, t121, t121; 0, -t150, 0, t150, -t147, -t147; t124, -t126, 0, t126, 0, 0; -t127, -t125, 0, t125, 0, 0; 0, -t129, 0, t129, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end