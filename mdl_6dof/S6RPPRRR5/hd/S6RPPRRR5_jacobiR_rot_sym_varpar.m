% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t5, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(4));
	t11 = sin(qJ(1));
	t16 = t11 * t10;
	t12 = cos(qJ(4));
	t15 = t11 * t12;
	t13 = cos(qJ(1));
	t14 = t13 * t10;
	t9 = t13 * t12;
	t1 = [-t16, 0, 0, t9, 0, 0; t14, 0, 0, t15, 0, 0; 0, 0, 0, -t10, 0, 0; -t15, 0, 0, -t14, 0, 0; t9, 0, 0, -t16, 0, 0; 0, 0, 0, -t12, 0, 0; -t13, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(4) + qJ(5);
	t22 = sin(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t22;
	t26 = cos(qJ(1));
	t27 = t26 * t22;
	t23 = cos(t24);
	t21 = t26 * t23;
	t20 = t25 * t23;
	t1 = [-t28, 0, 0, t21, t21, 0; t27, 0, 0, t20, t20, 0; 0, 0, 0, -t22, -t22, 0; -t20, 0, 0, -t27, -t27, 0; t21, 0, 0, -t28, -t28, 0; 0, 0, 0, -t23, -t23, 0; -t26, 0, 0, 0, 0, 0; -t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t105 = qJ(4) + qJ(5);
	t103 = sin(t105);
	t108 = cos(qJ(6));
	t116 = t103 * t108;
	t106 = sin(qJ(6));
	t107 = sin(qJ(1));
	t115 = t107 * t106;
	t114 = t107 * t108;
	t109 = cos(qJ(1));
	t113 = t109 * t106;
	t112 = t109 * t108;
	t104 = cos(t105);
	t111 = t104 * t115;
	t110 = t104 * t113;
	t102 = t109 * t103;
	t101 = t107 * t103;
	t100 = t103 * t106;
	t99 = t104 * t112;
	t98 = t104 * t114;
	t97 = t103 * t112 - t115;
	t96 = -t103 * t113 - t114;
	t95 = -t103 * t114 - t113;
	t94 = t103 * t115 - t112;
	t1 = [t95, 0, 0, t99, t99, t96; t97, 0, 0, t98, t98, -t94; 0, 0, 0, -t116, -t116, -t104 * t106; t94, 0, 0, -t110, -t110, -t97; t96, 0, 0, -t111, -t111, t95; 0, 0, 0, t100, t100, -t104 * t108; t107 * t104, 0, 0, t102, t102, 0; -t109 * t104, 0, 0, t101, t101, 0; 0, 0, 0, t104, t104, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end