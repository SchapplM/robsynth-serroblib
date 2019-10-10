% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t30 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t30, 0, 0, 0, 0; 0, -cos(pkin(11)) * t30, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t71 = sin(pkin(12));
	t74 = cos(pkin(12));
	t77 = sin(qJ(2));
	t78 = cos(qJ(2));
	t79 = t71 * t77 - t74 * t78;
	t76 = cos(pkin(6));
	t75 = cos(pkin(11));
	t73 = sin(pkin(6));
	t72 = sin(pkin(11));
	t70 = -t78 * t71 - t77 * t74;
	t69 = t79 * t76;
	t1 = [0, t72 * t73, 0, -t72 * t69 - t75 * t70, 0, 0; 0, -t75 * t73, 0, t75 * t69 - t72 * t70, 0, 0; 0, t76, 0, t79 * t73, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->5), mult. (42->12), div. (0->0), fcn. (65->8), ass. (0->15)
	t90 = sin(pkin(12));
	t93 = cos(pkin(12));
	t96 = sin(qJ(2));
	t97 = cos(qJ(2));
	t98 = t90 * t96 - t93 * t97;
	t95 = cos(pkin(6));
	t94 = cos(pkin(11));
	t92 = sin(pkin(6));
	t91 = sin(pkin(11));
	t89 = -t97 * t90 - t96 * t93;
	t88 = t98 * t95;
	t87 = t98 * t92;
	t86 = -t91 * t88 - t94 * t89;
	t85 = t94 * t88 - t91 * t89;
	t1 = [0, t91 * t92, 0, t86, t86, 0; 0, -t94 * t92, 0, t85, t85, 0; 0, t95, 0, t87, t87, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->11), mult. (70->24), div. (0->0), fcn. (106->10), ass. (0->21)
	t143 = sin(pkin(11));
	t144 = sin(pkin(6));
	t153 = t143 * t144;
	t146 = cos(pkin(11));
	t152 = t146 * t144;
	t142 = sin(pkin(12));
	t145 = cos(pkin(12));
	t148 = sin(qJ(2));
	t149 = cos(qJ(2));
	t151 = t149 * t142 + t148 * t145;
	t150 = t148 * t142 - t149 * t145;
	t147 = cos(pkin(6));
	t141 = qJ(4) + qJ(5);
	t140 = cos(t141);
	t139 = sin(t141);
	t136 = t151 * t147;
	t135 = t150 * t147;
	t134 = t150 * t144;
	t133 = -t143 * t135 + t146 * t151;
	t132 = t146 * t135 + t143 * t151;
	t1 = [0, t153, 0, t133, t133, (-t143 * t136 - t146 * t150) * t139 - t140 * t153; 0, -t152, 0, t132, t132, (t146 * t136 - t143 * t150) * t139 + t140 * t152; 0, t147, 0, t134, t134, t151 * t139 * t144 - t147 * t140;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end