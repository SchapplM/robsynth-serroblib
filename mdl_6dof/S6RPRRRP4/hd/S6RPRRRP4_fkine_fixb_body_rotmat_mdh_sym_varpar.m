% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t63 = cos(pkin(10));
	t62 = sin(pkin(10));
	t1 = [t65 * t63, -t65 * t62, t64, t65 * pkin(1) + t64 * qJ(2) + 0; t64 * t63, -t64 * t62, -t65, t64 * pkin(1) - t65 * qJ(2) + 0; t62, t63, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = pkin(7) + qJ(2);
	t69 = pkin(10) + qJ(3);
	t68 = cos(t69);
	t67 = sin(t69);
	t66 = cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t72 * t68, -t72 * t67, t71, t72 * t66 + t70 * t71 + 0; t71 * t68, -t71 * t67, -t72, t71 * t66 - t72 * t70 + 0; t67, t68, 0, sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t78 = pkin(10) + qJ(3);
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t77 = -pkin(8) - pkin(7) - qJ(2);
	t76 = qJ(4) + t78;
	t75 = cos(t76);
	t74 = sin(t76);
	t73 = pkin(3) * cos(t78) + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t80 * t75, -t80 * t74, t79, t80 * t73 - t79 * t77 + 0; t79 * t75, -t79 * t74, -t80, t79 * t73 + t80 * t77 + 0; t74, t75, 0, pkin(3) * sin(t78) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t87 = sin(qJ(5));
	t88 = sin(qJ(1));
	t95 = t88 * t87;
	t89 = cos(qJ(5));
	t94 = t88 * t89;
	t90 = cos(qJ(1));
	t93 = t90 * t87;
	t92 = t90 * t89;
	t86 = pkin(10) + qJ(3);
	t84 = qJ(4) + t86;
	t82 = sin(t84);
	t83 = cos(t84);
	t91 = pkin(4) * t83 + pkin(9) * t82 + pkin(3) * cos(t86) + cos(pkin(10)) * pkin(2) + pkin(1);
	t85 = -pkin(8) - pkin(7) - qJ(2);
	t1 = [t83 * t92 + t95, -t83 * t93 + t94, t90 * t82, -t88 * t85 + t91 * t90 + 0; t83 * t94 - t93, -t83 * t95 - t92, t88 * t82, t90 * t85 + t91 * t88 + 0; t82 * t89, -t82 * t87, -t83, t82 * pkin(4) - t83 * pkin(9) + pkin(3) * sin(t86) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:51:58
	% EndTime: 2020-11-04 21:51:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->26), mult. (43->26), div. (0->0), fcn. (56->10), ass. (0->17)
	t105 = sin(qJ(5));
	t106 = sin(qJ(1));
	t113 = t106 * t105;
	t107 = cos(qJ(5));
	t112 = t106 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t105;
	t110 = t108 * t107;
	t103 = pkin(10) + qJ(3);
	t100 = t107 * pkin(5) + pkin(4);
	t104 = qJ(6) + pkin(9);
	t102 = qJ(4) + t103;
	t97 = sin(t102);
	t98 = cos(t102);
	t109 = pkin(3) * cos(t103) + t100 * t98 + t104 * t97 + cos(pkin(10)) * pkin(2) + pkin(1);
	t96 = t105 * pkin(5) + pkin(7) + pkin(8) + qJ(2);
	t1 = [t98 * t110 + t113, -t98 * t111 + t112, t108 * t97, t96 * t106 + t109 * t108 + 0; t98 * t112 - t111, -t98 * t113 - t110, t106 * t97, t109 * t106 - t96 * t108 + 0; t97 * t107, -t97 * t105, -t98, t97 * t100 - t98 * t104 + pkin(3) * sin(t103) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end