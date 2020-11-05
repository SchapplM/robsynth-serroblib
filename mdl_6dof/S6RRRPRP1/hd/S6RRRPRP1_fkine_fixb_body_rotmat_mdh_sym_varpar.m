% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:26
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
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
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
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
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t73 = qJ(2) + qJ(3);
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t72 = -qJ(4) - pkin(8) - pkin(7);
	t71 = pkin(10) + t73;
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = pkin(3) * cos(t73) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t70, -t75 * t69, t74, t75 * t68 - t74 * t72 + 0; t74 * t70, -t74 * t69, -t75, t74 * t68 + t75 * t72 + 0; t69, t70, 0, pkin(3) * sin(t73) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t82 = sin(qJ(5));
	t83 = sin(qJ(1));
	t90 = t83 * t82;
	t84 = cos(qJ(5));
	t89 = t83 * t84;
	t85 = cos(qJ(1));
	t88 = t85 * t82;
	t87 = t85 * t84;
	t81 = qJ(2) + qJ(3);
	t79 = pkin(10) + t81;
	t77 = sin(t79);
	t78 = cos(t79);
	t86 = pkin(4) * t78 + pkin(9) * t77 + pkin(3) * cos(t81) + cos(qJ(2)) * pkin(2) + pkin(1);
	t80 = -qJ(4) - pkin(8) - pkin(7);
	t1 = [t78 * t87 + t90, -t78 * t88 + t89, t85 * t77, -t83 * t80 + t86 * t85 + 0; t78 * t89 - t88, -t78 * t90 - t87, t83 * t77, t85 * t80 + t86 * t83 + 0; t77 * t84, -t77 * t82, -t78, t77 * pkin(4) - t78 * pkin(9) + pkin(3) * sin(t81) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:26:02
	% EndTime: 2020-11-04 22:26:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->26), mult. (43->26), div. (0->0), fcn. (56->10), ass. (0->17)
	t100 = sin(qJ(1));
	t99 = sin(qJ(5));
	t108 = t100 * t99;
	t102 = cos(qJ(1));
	t107 = t102 * t99;
	t101 = cos(qJ(5));
	t106 = t100 * t101;
	t105 = t102 * t101;
	t97 = qJ(2) + qJ(3);
	t104 = pkin(5) * t99 + pkin(7) + pkin(8) + qJ(4);
	t95 = pkin(10) + t97;
	t92 = sin(t95);
	t93 = cos(t95);
	t94 = t101 * pkin(5) + pkin(4);
	t98 = -qJ(6) - pkin(9);
	t103 = -t92 * t98 + t93 * t94 + pkin(3) * cos(t97) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t93 * t105 + t108, -t93 * t107 + t106, t102 * t92, t104 * t100 + t103 * t102 + 0; t93 * t106 - t107, -t93 * t108 - t105, t100 * t92, t103 * t100 - t104 * t102 + 0; t92 * t101, -t92 * t99, -t93, t92 * t94 + t93 * t98 + pkin(3) * sin(t97) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end