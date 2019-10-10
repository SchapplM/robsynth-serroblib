% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t60 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t60, 0, 0, 0, 0; 0, -cos(qJ(1)) * t60, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t74 = cos(pkin(6));
	t77 = cos(qJ(2));
	t79 = t74 * t77;
	t78 = cos(qJ(1));
	t76 = sin(qJ(1));
	t75 = sin(qJ(2));
	t73 = sin(pkin(6));
	t1 = [0, t76 * t73, 0, t78 * t75 + t76 * t79, 0, 0; 0, -t78 * t73, 0, t76 * t75 - t78 * t79, 0, 0; 1, t74, 0, -t73 * t77, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t111 = sin(pkin(6));
	t114 = sin(qJ(1));
	t122 = t114 * t111;
	t113 = sin(qJ(2));
	t121 = t114 * t113;
	t115 = cos(qJ(2));
	t120 = t114 * t115;
	t116 = cos(qJ(1));
	t119 = t116 * t111;
	t118 = t116 * t113;
	t117 = t116 * t115;
	t112 = cos(pkin(6));
	t110 = pkin(11) + qJ(4);
	t109 = cos(t110);
	t108 = sin(t110);
	t1 = [0, t122, 0, t112 * t120 + t118, (-t112 * t121 + t117) * t108 - t109 * t122, 0; 0, -t119, 0, -t112 * t117 + t121, (t112 * t118 + t120) * t108 + t109 * t119, 0; 1, t112, 0, -t111 * t115, t108 * t111 * t113 - t109 * t112, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t143 = sin(pkin(6));
	t146 = sin(qJ(1));
	t154 = t146 * t143;
	t145 = sin(qJ(2));
	t153 = t146 * t145;
	t147 = cos(qJ(2));
	t152 = t146 * t147;
	t148 = cos(qJ(1));
	t151 = t148 * t143;
	t150 = t148 * t145;
	t149 = t148 * t147;
	t144 = cos(pkin(6));
	t142 = pkin(11) + qJ(4);
	t141 = cos(t142);
	t140 = sin(t142);
	t1 = [0, t154, 0, t144 * t152 + t150, (-t144 * t153 + t149) * t140 - t141 * t154, 0; 0, -t151, 0, -t144 * t149 + t153, (t144 * t150 + t152) * t140 + t141 * t151, 0; 1, t144, 0, -t143 * t147, t143 * t145 * t140 - t144 * t141, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end