% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, -cos(qJ(1)) * t59, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t91 = sin(pkin(12));
	t93 = cos(pkin(12));
	t95 = sin(qJ(2));
	t97 = cos(qJ(2));
	t99 = t91 * t95 - t93 * t97;
	t98 = cos(qJ(1));
	t96 = sin(qJ(1));
	t94 = cos(pkin(6));
	t92 = sin(pkin(6));
	t90 = -t97 * t91 - t95 * t93;
	t89 = t99 * t94;
	t1 = [0, t96 * t92, 0, -t96 * t89 - t98 * t90, 0, 0; 0, -t98 * t92, 0, t98 * t89 - t96 * t90, 0, 0; 1, t94, 0, t99 * t92, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (15->5), mult. (42->12), div. (0->0), fcn. (65->8), ass. (0->15)
	t109 = sin(pkin(12));
	t111 = cos(pkin(12));
	t113 = sin(qJ(2));
	t115 = cos(qJ(2));
	t117 = t109 * t113 - t111 * t115;
	t116 = cos(qJ(1));
	t114 = sin(qJ(1));
	t112 = cos(pkin(6));
	t110 = sin(pkin(6));
	t108 = -t115 * t109 - t113 * t111;
	t107 = t117 * t112;
	t106 = t117 * t110;
	t105 = -t114 * t107 - t116 * t108;
	t104 = t116 * t107 - t114 * t108;
	t1 = [0, t114 * t110, 0, t105, t105, 0; 0, -t116 * t110, 0, t104, t104, 0; 1, t112, 0, t106, t106, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->11), mult. (70->24), div. (0->0), fcn. (106->10), ass. (0->21)
	t159 = sin(pkin(6));
	t163 = sin(qJ(1));
	t169 = t163 * t159;
	t165 = cos(qJ(1));
	t168 = t165 * t159;
	t158 = sin(pkin(12));
	t160 = cos(pkin(12));
	t162 = sin(qJ(2));
	t164 = cos(qJ(2));
	t167 = t164 * t158 + t162 * t160;
	t166 = t162 * t158 - t164 * t160;
	t161 = cos(pkin(6));
	t157 = qJ(4) + qJ(5);
	t156 = cos(t157);
	t155 = sin(t157);
	t152 = t167 * t161;
	t151 = t166 * t161;
	t150 = t166 * t159;
	t149 = -t163 * t151 + t165 * t167;
	t148 = t165 * t151 + t163 * t167;
	t1 = [0, t169, 0, t149, t149, (-t163 * t152 - t165 * t166) * t155 - t156 * t169; 0, -t168, 0, t148, t148, (t165 * t152 - t163 * t166) * t155 + t156 * t168; 1, t161, 0, t150, t150, t167 * t155 * t159 - t161 * t156;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end