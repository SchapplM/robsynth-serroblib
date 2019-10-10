% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPPRR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t30 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t30, 0, 0, 0, 0; 0, -cos(pkin(10)) * t30, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t54 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t54, 0, 0, 0, 0; 0, -cos(pkin(10)) * t54, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t74 = sin(pkin(11));
	t77 = cos(pkin(11));
	t80 = sin(qJ(2));
	t81 = cos(qJ(2));
	t82 = t74 * t81 + t77 * t80;
	t79 = cos(pkin(6));
	t78 = cos(pkin(10));
	t76 = sin(pkin(6));
	t75 = sin(pkin(10));
	t73 = -t80 * t74 + t81 * t77;
	t72 = t82 * t79;
	t1 = [0, t75 * t76, 0, 0, -t75 * t72 + t78 * t73, 0; 0, -t78 * t76, 0, 0, t78 * t72 + t75 * t73, 0; 0, t79, 0, 0, t82 * t76, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t114 = sin(pkin(10));
	t115 = sin(pkin(6));
	t126 = t114 * t115;
	t117 = cos(pkin(10));
	t125 = t117 * t115;
	t113 = sin(pkin(11));
	t116 = cos(pkin(11));
	t120 = sin(qJ(2));
	t122 = cos(qJ(2));
	t124 = t122 * t113 + t120 * t116;
	t123 = t120 * t113 - t122 * t116;
	t121 = cos(qJ(5));
	t119 = sin(qJ(5));
	t118 = cos(pkin(6));
	t110 = t124 * t118;
	t109 = t123 * t118;
	t1 = [0, t126, 0, 0, -t114 * t110 - t117 * t123, t119 * t126 - (-t114 * t109 + t117 * t124) * t121; 0, -t125, 0, 0, t117 * t110 - t114 * t123, -t119 * t125 - (t117 * t109 + t114 * t124) * t121; 0, t118, 0, 0, t124 * t115, -t123 * t121 * t115 + t118 * t119;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end