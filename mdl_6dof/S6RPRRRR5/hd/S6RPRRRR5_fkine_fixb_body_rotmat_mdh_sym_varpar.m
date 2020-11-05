% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = cos(pkin(11));
	t65 = sin(pkin(11));
	t1 = [t68 * t66, -t68 * t65, t67, t68 * pkin(1) + t67 * qJ(2) + 0; t67 * t66, -t67 * t65, -t68, t67 * pkin(1) - t68 * qJ(2) + 0; t65, t66, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = pkin(7) + qJ(2);
	t72 = pkin(11) + qJ(3);
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = cos(pkin(11)) * pkin(2) + pkin(1);
	t1 = [t75 * t71, -t75 * t70, t74, t75 * t69 + t73 * t74 + 0; t74 * t71, -t74 * t70, -t75, t74 * t69 - t75 * t73 + 0; t70, t71, 0, sin(pkin(11)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t81 = pkin(11) + qJ(3);
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t80 = -pkin(8) - pkin(7) - qJ(2);
	t79 = qJ(4) + t81;
	t78 = cos(t79);
	t77 = sin(t79);
	t76 = pkin(3) * cos(t81) + cos(pkin(11)) * pkin(2) + pkin(1);
	t1 = [t83 * t78, -t83 * t77, t82, t83 * t76 - t82 * t80 + 0; t82 * t78, -t82 * t77, -t83, t82 * t76 + t83 * t80 + 0; t77, t78, 0, pkin(3) * sin(t81) + sin(pkin(11)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t90 = sin(qJ(5));
	t91 = sin(qJ(1));
	t98 = t91 * t90;
	t92 = cos(qJ(5));
	t97 = t91 * t92;
	t93 = cos(qJ(1));
	t96 = t93 * t90;
	t95 = t93 * t92;
	t89 = pkin(11) + qJ(3);
	t87 = qJ(4) + t89;
	t85 = sin(t87);
	t86 = cos(t87);
	t94 = pkin(4) * t86 + pkin(9) * t85 + pkin(3) * cos(t89) + cos(pkin(11)) * pkin(2) + pkin(1);
	t88 = -pkin(8) - pkin(7) - qJ(2);
	t1 = [t86 * t95 + t98, -t86 * t96 + t97, t93 * t85, -t91 * t88 + t94 * t93 + 0; t86 * t97 - t96, -t86 * t98 - t95, t91 * t85, t93 * t88 + t94 * t91 + 0; t85 * t92, -t85 * t90, -t86, t85 * pkin(4) - t86 * pkin(9) + pkin(3) * sin(t89) + sin(pkin(11)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:11
	% EndTime: 2020-11-04 21:56:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->27), mult. (43->26), div. (0->0), fcn. (56->12), ass. (0->18)
	t109 = qJ(5) + qJ(6);
	t106 = sin(t109);
	t110 = sin(qJ(1));
	t117 = t110 * t106;
	t107 = cos(t109);
	t116 = t110 * t107;
	t111 = cos(qJ(1));
	t115 = t111 * t106;
	t114 = t111 * t107;
	t108 = pkin(11) + qJ(3);
	t105 = qJ(4) + t108;
	t100 = sin(t105);
	t101 = cos(t105);
	t103 = cos(qJ(5)) * pkin(5) + pkin(4);
	t112 = pkin(10) + pkin(9);
	t113 = pkin(3) * cos(t108) + t100 * t112 + t101 * t103 + cos(pkin(11)) * pkin(2) + pkin(1);
	t99 = sin(qJ(5)) * pkin(5) + qJ(2) + pkin(7) + pkin(8);
	t1 = [t101 * t114 + t117, -t101 * t115 + t116, t111 * t100, t99 * t110 + t113 * t111 + 0; t101 * t116 - t115, -t101 * t117 - t114, t110 * t100, t113 * t110 - t99 * t111 + 0; t100 * t107, -t100 * t106, -t101, t100 * t103 - t101 * t112 + pkin(3) * sin(t108) + sin(pkin(11)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end