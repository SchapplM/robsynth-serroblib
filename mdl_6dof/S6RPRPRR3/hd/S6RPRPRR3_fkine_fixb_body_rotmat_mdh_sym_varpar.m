% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t61 = qJ(1) + pkin(10);
	t60 = cos(t61);
	t59 = sin(t61);
	t1 = [t60, -t59, 0, cos(qJ(1)) * pkin(1) + 0; t59, t60, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t66 = cos(qJ(3));
	t65 = sin(qJ(3));
	t64 = qJ(1) + pkin(10);
	t63 = cos(t64);
	t62 = sin(t64);
	t1 = [t63 * t66, -t63 * t65, t62, t63 * pkin(2) + t62 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t62 * t66, -t62 * t65, -t63, t62 * pkin(2) - t63 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t65, t66, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t70 = sin(pkin(11));
	t73 = cos(qJ(3));
	t76 = t70 * t73;
	t71 = cos(pkin(11));
	t75 = t71 * t73;
	t72 = sin(qJ(3));
	t74 = pkin(3) * t73 + qJ(4) * t72 + pkin(2);
	t69 = qJ(1) + pkin(10);
	t68 = cos(t69);
	t67 = sin(t69);
	t1 = [t67 * t70 + t68 * t75, t67 * t71 - t68 * t76, t68 * t72, cos(qJ(1)) * pkin(1) + t67 * pkin(7) + 0 + t74 * t68; t67 * t75 - t68 * t70, -t67 * t76 - t68 * t71, t67 * t72, sin(qJ(1)) * pkin(1) - t68 * pkin(7) + 0 + t74 * t67; t72 * t71, -t72 * t70, -t73, t72 * pkin(3) - t73 * qJ(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t83 = qJ(1) + pkin(10);
	t79 = sin(t83);
	t87 = cos(qJ(3));
	t91 = t79 * t87;
	t81 = cos(t83);
	t90 = t81 * t87;
	t89 = sin(pkin(11)) * pkin(4) + pkin(7);
	t77 = cos(pkin(11)) * pkin(4) + pkin(3);
	t85 = -pkin(8) - qJ(4);
	t86 = sin(qJ(3));
	t88 = t77 * t87 - t85 * t86 + pkin(2);
	t82 = pkin(11) + qJ(5);
	t80 = cos(t82);
	t78 = sin(t82);
	t1 = [t79 * t78 + t80 * t90, -t78 * t90 + t79 * t80, t81 * t86, cos(qJ(1)) * pkin(1) + 0 + t89 * t79 + t88 * t81; -t81 * t78 + t80 * t91, -t78 * t91 - t81 * t80, t79 * t86, sin(qJ(1)) * pkin(1) + 0 - t89 * t81 + t88 * t79; t86 * t80, -t86 * t78, -t87, t86 * t77 + t87 * t85 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:39:50
	% EndTime: 2020-11-04 21:39:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->26), mult. (44->28), div. (0->0), fcn. (57->12), ass. (0->16)
	t100 = pkin(11) + qJ(5);
	t107 = pkin(7) + pkin(5) * sin(t100) + sin(pkin(11)) * pkin(4);
	t103 = cos(qJ(3));
	t101 = qJ(1) + pkin(10);
	t96 = sin(t101);
	t106 = t103 * t96;
	t97 = cos(t101);
	t105 = t103 * t97;
	t102 = sin(qJ(3));
	t92 = pkin(5) * cos(t100) + cos(pkin(11)) * pkin(4) + pkin(3);
	t99 = -pkin(9) - pkin(8) - qJ(4);
	t104 = -t102 * t99 + t103 * t92 + pkin(2);
	t98 = qJ(6) + t100;
	t95 = cos(t98);
	t94 = sin(t98);
	t1 = [t95 * t105 + t96 * t94, -t94 * t105 + t96 * t95, t97 * t102, cos(qJ(1)) * pkin(1) + 0 + t107 * t96 + t104 * t97; t95 * t106 - t97 * t94, -t94 * t106 - t97 * t95, t96 * t102, sin(qJ(1)) * pkin(1) + 0 - t107 * t97 + t104 * t96; t102 * t95, -t102 * t94, -t103, t102 * t92 + t103 * t99 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end