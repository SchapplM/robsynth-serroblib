% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(11));
	t55 = sin(pkin(6));
	t54 = sin(pkin(11));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0, 0; 0, t57, -t55 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t64 = cos(pkin(6));
	t66 = cos(qJ(2));
	t67 = t64 * t66;
	t65 = sin(qJ(2));
	t63 = cos(pkin(11));
	t62 = sin(pkin(6));
	t61 = sin(pkin(11));
	t1 = [0, t61 * t62, t61 * t67 + t63 * t65, 0, 0, 0; 0, -t63 * t62, t61 * t65 - t63 * t67, 0, 0, 0; 0, t64, -t62 * t66, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t76 = sin(pkin(6));
	t80 = cos(qJ(2));
	t82 = t76 * t80;
	t78 = cos(pkin(6));
	t81 = t78 * t80;
	t79 = sin(qJ(2));
	t77 = cos(pkin(11));
	t75 = sin(pkin(11));
	t74 = t75 * t81 + t77 * t79;
	t73 = t75 * t79 - t77 * t81;
	t1 = [0, t75 * t76, t74, 0, t74, 0; 0, -t77 * t76, t73, 0, t73, 0; 0, t78, -t82, 0, -t82, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
	t121 = sin(pkin(11));
	t122 = sin(pkin(6));
	t131 = t121 * t122;
	t126 = cos(qJ(2));
	t130 = t122 * t126;
	t123 = cos(pkin(11));
	t129 = t123 * t122;
	t124 = cos(pkin(6));
	t125 = sin(qJ(2));
	t128 = t124 * t125;
	t127 = t124 * t126;
	t120 = qJ(3) + pkin(12) + qJ(5);
	t119 = cos(t120);
	t118 = sin(t120);
	t117 = t121 * t127 + t123 * t125;
	t116 = t121 * t125 - t123 * t127;
	t1 = [0, t131, t117, 0, t117, (-t121 * t128 + t123 * t126) * t118 - t119 * t131; 0, -t129, t116, 0, t116, (t121 * t126 + t123 * t128) * t118 + t119 * t129; 0, t124, -t130, 0, -t130, t122 * t125 * t118 - t124 * t119;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end