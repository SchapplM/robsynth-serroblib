% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:49
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRP10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t70 = cos(pkin(6));
	t73 = cos(qJ(2));
	t75 = t70 * t73;
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t69 = sin(pkin(6));
	t1 = [0, t72 * t69, t74 * t71 + t72 * t75, 0, 0, 0; 0, -t74 * t69, t72 * t71 - t74 * t75, 0, 0, 0; 1, t70, -t69 * t73, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:38
	% EndTime: 2019-10-10 11:49:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t95 = cos(pkin(6));
	t98 = cos(qJ(2));
	t100 = t95 * t98;
	t99 = cos(qJ(1));
	t97 = sin(qJ(1));
	t96 = sin(qJ(2));
	t94 = sin(pkin(6));
	t1 = [0, t97 * t94, t97 * t100 + t99 * t96, 0, 0, 0; 0, -t99 * t94, -t99 * t100 + t97 * t96, 0, 0, 0; 1, t95, -t94 * t98, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:38
	% EndTime: 2019-10-10 11:49:38
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t103 = sin(pkin(6));
	t107 = sin(qJ(1));
	t116 = t107 * t103;
	t106 = sin(qJ(2));
	t115 = t107 * t106;
	t109 = cos(qJ(2));
	t114 = t107 * t109;
	t110 = cos(qJ(1));
	t113 = t110 * t103;
	t112 = t110 * t106;
	t111 = t110 * t109;
	t108 = cos(qJ(3));
	t105 = sin(qJ(3));
	t104 = cos(pkin(6));
	t1 = [0, t116, t104 * t114 + t112, 0, (-t104 * t115 + t111) * t105 - t108 * t116, 0; 0, -t113, -t104 * t111 + t115, 0, (t104 * t112 + t114) * t105 + t108 * t113, 0; 1, t104, -t103 * t109, 0, t103 * t106 * t105 - t104 * t108, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:38
	% EndTime: 2019-10-10 11:49:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t140 = sin(pkin(6));
	t144 = sin(qJ(1));
	t153 = t144 * t140;
	t143 = sin(qJ(2));
	t152 = t144 * t143;
	t146 = cos(qJ(2));
	t151 = t144 * t146;
	t147 = cos(qJ(1));
	t150 = t147 * t140;
	t149 = t147 * t143;
	t148 = t147 * t146;
	t145 = cos(qJ(3));
	t142 = sin(qJ(3));
	t141 = cos(pkin(6));
	t1 = [0, t153, t141 * t151 + t149, 0, (-t141 * t152 + t148) * t142 - t145 * t153, 0; 0, -t150, -t141 * t148 + t152, 0, (t141 * t149 + t151) * t142 + t145 * t150, 0; 1, t141, -t140 * t146, 0, t140 * t143 * t142 - t141 * t145, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end