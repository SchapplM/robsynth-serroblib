% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t1 = [t60 * t59, -t60 * t57, t58, t60 * pkin(1) + t58 * pkin(7) + 0; t58 * t59, -t58 * t57, -t60, t58 * pkin(1) - t60 * pkin(7) + 0; t57, t59, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t67 = pkin(8) + pkin(7);
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t64 = qJ(2) + qJ(3);
	t63 = cos(t64);
	t62 = sin(t64);
	t61 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t66 * t63, -t66 * t62, t65, t66 * t61 + t67 * t65 + 0; t65 * t63, -t65 * t62, -t66, t65 * t61 - t66 * t67 + 0; t62, t63, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t73 = qJ(2) + qJ(3);
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t72 = -qJ(4) - pkin(8) - pkin(7);
	t71 = pkin(11) + t73;
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = pkin(3) * cos(t73) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t70, -t75 * t69, t74, t75 * t68 - t74 * t72 + 0; t74 * t70, -t74 * t69, -t75, t74 * t68 + t75 * t72 + 0; t69, t70, 0, pkin(3) * sin(t73) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (50->18), mult. (17->14), div. (0->0), fcn. (25->10), ass. (0->10)
	t82 = qJ(2) + qJ(3);
	t80 = pkin(11) + t82;
	t84 = cos(qJ(1));
	t83 = sin(qJ(1));
	t81 = -pkin(9) - qJ(4) - pkin(8) - pkin(7);
	t79 = qJ(5) + t80;
	t78 = cos(t79);
	t77 = sin(t79);
	t76 = pkin(4) * cos(t80) + pkin(3) * cos(t82) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t84 * t78, -t84 * t77, t83, t84 * t76 - t83 * t81 + 0; t83 * t78, -t83 * t77, -t84, t83 * t76 + t84 * t81 + 0; t77, t78, 0, pkin(4) * sin(t80) + pkin(3) * sin(t82) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:59
	% EndTime: 2020-11-04 22:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (86->25), mult. (39->26), div. (0->0), fcn. (52->12), ass. (0->16)
	t92 = sin(qJ(6));
	t93 = sin(qJ(1));
	t100 = t93 * t92;
	t94 = cos(qJ(6));
	t99 = t93 * t94;
	t95 = cos(qJ(1));
	t98 = t95 * t92;
	t97 = t95 * t94;
	t91 = qJ(2) + qJ(3);
	t89 = pkin(11) + t91;
	t88 = qJ(5) + t89;
	t86 = sin(t88);
	t87 = cos(t88);
	t96 = pkin(5) * t87 + pkin(10) * t86 + pkin(4) * cos(t89) + pkin(3) * cos(t91) + cos(qJ(2)) * pkin(2) + pkin(1);
	t90 = -pkin(9) - qJ(4) - pkin(8) - pkin(7);
	t1 = [t87 * t97 + t100, -t87 * t98 + t99, t95 * t86, -t93 * t90 + t96 * t95 + 0; t87 * t99 - t98, -t87 * t100 - t97, t93 * t86, t95 * t90 + t96 * t93 + 0; t86 * t94, -t86 * t92, -t87, t86 * pkin(5) - t87 * pkin(10) + pkin(4) * sin(t89) + pkin(3) * sin(t91) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end