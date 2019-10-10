% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t30 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t30, 0, 0, 0, 0; 0, -cos(pkin(10)) * t30, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t71 = sin(pkin(11));
	t74 = cos(pkin(11));
	t77 = sin(qJ(2));
	t78 = cos(qJ(2));
	t79 = t71 * t77 - t74 * t78;
	t76 = cos(pkin(6));
	t75 = cos(pkin(10));
	t73 = sin(pkin(6));
	t72 = sin(pkin(10));
	t70 = -t78 * t71 - t77 * t74;
	t69 = t79 * t76;
	t1 = [0, t72 * t73, 0, -t72 * t69 - t75 * t70, 0, 0; 0, -t75 * t73, 0, t75 * t69 - t72 * t70, 0, 0; 0, t76, 0, t79 * t73, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t111 = sin(pkin(10));
	t112 = sin(pkin(6));
	t123 = t111 * t112;
	t114 = cos(pkin(10));
	t122 = t114 * t112;
	t110 = sin(pkin(11));
	t113 = cos(pkin(11));
	t117 = sin(qJ(2));
	t119 = cos(qJ(2));
	t121 = t110 * t119 + t113 * t117;
	t120 = t110 * t117 - t113 * t119;
	t118 = cos(qJ(4));
	t116 = sin(qJ(4));
	t115 = cos(pkin(6));
	t107 = t121 * t115;
	t106 = t120 * t115;
	t1 = [0, t123, 0, -t106 * t111 + t114 * t121, (-t107 * t111 - t114 * t120) * t116 - t118 * t123, 0; 0, -t122, 0, t106 * t114 + t111 * t121, (t107 * t114 - t111 * t120) * t116 + t118 * t122, 0; 0, t115, 0, t120 * t112, t112 * t116 * t121 - t115 * t118, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t136 = sin(pkin(10));
	t137 = sin(pkin(6));
	t148 = t136 * t137;
	t139 = cos(pkin(10));
	t147 = t139 * t137;
	t135 = sin(pkin(11));
	t138 = cos(pkin(11));
	t142 = sin(qJ(2));
	t144 = cos(qJ(2));
	t146 = t144 * t135 + t142 * t138;
	t145 = t142 * t135 - t144 * t138;
	t143 = cos(qJ(4));
	t141 = sin(qJ(4));
	t140 = cos(pkin(6));
	t132 = t146 * t140;
	t131 = t145 * t140;
	t1 = [0, t148, 0, -t136 * t131 + t139 * t146, (-t136 * t132 - t139 * t145) * t141 - t143 * t148, 0; 0, -t147, 0, t139 * t131 + t136 * t146, (t139 * t132 - t136 * t145) * t141 + t143 * t147, 0; 0, t140, 0, t145 * t137, t146 * t141 * t137 - t140 * t143, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end