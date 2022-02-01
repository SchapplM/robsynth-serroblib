% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   Siehe auch: S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(7));
	t7 = sin(pkin(7));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(8));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(8));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(7));
	t39 = sin(pkin(7));
	t1 = [-t41 * t46 + t45, 0, 0, 0, 0; t41 * t44 + t47, 0, 0, 0, 0; 0, 0, 0, 0, 0; t41 * t47 + t44, 0, 0, 0, 0; -t41 * t45 + t46, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t39 * t42, 0, 0, 0, 0; t39 * t43, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (30->18), div. (0->0), fcn. (46->8), ass. (0->17)
	t61 = sin(pkin(7));
	t65 = sin(qJ(1));
	t72 = t61 * t65;
	t66 = cos(qJ(1));
	t71 = t61 * t66;
	t60 = sin(pkin(8));
	t70 = t65 * t60;
	t63 = cos(pkin(8));
	t69 = t65 * t63;
	t68 = t66 * t60;
	t67 = t66 * t63;
	t64 = cos(pkin(7));
	t62 = cos(pkin(9));
	t59 = sin(pkin(9));
	t58 = t64 * t67 + t70;
	t57 = -t64 * t69 + t68;
	t1 = [t57 * t62 - t59 * t72, 0, 0, 0, 0; t58 * t62 + t59 * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t57 * t59 - t62 * t72, 0, 0, 0, 0; -t58 * t59 + t62 * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t64 * t70 - t67, 0, 0, 0, 0; t64 * t68 - t69, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (40->22), mult. (112->48), div. (0->0), fcn. (164->10), ass. (0->27)
	t106 = cos(pkin(9));
	t107 = cos(pkin(8));
	t109 = sin(qJ(5));
	t104 = sin(pkin(8));
	t111 = cos(qJ(5));
	t121 = t104 * t111;
	t101 = t106 * t121 - t109 * t107;
	t110 = sin(qJ(1));
	t112 = cos(qJ(1));
	t108 = cos(pkin(7));
	t122 = t104 * t109;
	t103 = sin(pkin(9));
	t105 = sin(pkin(7));
	t119 = t106 * t107;
	t99 = t105 * t103 + t108 * t119;
	t113 = t108 * t122 + t99 * t111;
	t125 = t112 * t101 - t113 * t110;
	t100 = t106 * t122 + t111 * t107;
	t95 = t108 * t121 - t99 * t109;
	t124 = -t110 * t100 + t95 * t112;
	t120 = t105 * t106;
	t118 = t107 * t110;
	t116 = t110 * t104;
	t115 = t112 * t104;
	t114 = t112 * t107;
	t98 = -t108 * t103 + t105 * t119;
	t1 = [t125, 0, 0, 0, t124; (t106 * t116 + t99 * t112) * t111 + (t108 * t115 - t118) * t109, 0, 0, 0, -(-t106 * t115 + t99 * t110) * t109 + (t108 * t116 + t114) * t111; 0, 0, 0, 0, t105 * t121 - t98 * t109; -t112 * t100 - t95 * t110, 0, 0, 0, -t110 * t101 - t112 * t113; t124, 0, 0, 0, t125; 0, 0, 0, 0, -t105 * t122 - t98 * t111; (-t108 * t118 + t115) * t103 + t110 * t120, 0, 0, 0, 0; (t108 * t114 + t116) * t103 - t112 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end