% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiR_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(2));
	t3 = sin(qJ(2));
	t1 = [0, -t3, 0, 0, 0; 0, 0, 0, 0, 0; 0, t4, 0, 0, 0; 0, -t4, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t3, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t20 = sin(qJ(3));
	t21 = sin(qJ(2));
	t27 = t21 * t20;
	t22 = cos(qJ(3));
	t26 = t21 * t22;
	t23 = cos(qJ(2));
	t25 = t23 * t20;
	t24 = t23 * t22;
	t1 = [0, -t26, -t25, 0, 0; 0, 0, -t22, 0, 0; 0, t24, -t27, 0, 0; 0, t27, -t24, 0, 0; 0, 0, t20, 0, 0; 0, -t25, -t26, 0, 0; 0, t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, t21, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t34 = qJ(3) + qJ(4);
	t32 = sin(t34);
	t35 = sin(qJ(2));
	t40 = t35 * t32;
	t33 = cos(t34);
	t39 = t35 * t33;
	t36 = cos(qJ(2));
	t38 = t36 * t32;
	t37 = t36 * t33;
	t1 = [0, -t39, -t38, -t38, 0; 0, 0, -t33, -t33, 0; 0, t37, -t40, -t40, 0; 0, t40, -t37, -t37, 0; 0, 0, t32, t32, 0; 0, -t38, -t39, -t39, 0; 0, t36, 0, 0, 0; 0, 0, 0, 0, 0; 0, t35, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t93 = qJ(3) + qJ(4);
	t92 = cos(t93);
	t96 = cos(qJ(5));
	t104 = t92 * t96;
	t94 = sin(qJ(5));
	t95 = sin(qJ(2));
	t103 = t95 * t94;
	t102 = t95 * t96;
	t97 = cos(qJ(2));
	t101 = t97 * t94;
	t100 = t97 * t96;
	t91 = sin(t93);
	t99 = t91 * t102;
	t98 = t91 * t100;
	t90 = t97 * t92;
	t89 = t92 * t94;
	t88 = t95 * t92;
	t87 = t91 * t101;
	t86 = t91 * t103;
	t85 = t92 * t100 + t103;
	t84 = -t92 * t101 + t102;
	t83 = -t92 * t102 + t101;
	t82 = t92 * t103 + t100;
	t1 = [0, t83, -t98, -t98, t84; 0, 0, -t104, -t104, t91 * t94; 0, t85, -t99, -t99, -t82; 0, t82, t87, t87, -t85; 0, 0, t89, t89, t91 * t96; 0, t84, t86, t86, t83; 0, -t95 * t91, t90, t90, 0; 0, 0, -t91, -t91, 0; 0, t97 * t91, t88, t88, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end