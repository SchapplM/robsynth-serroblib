% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t62, -t61, 0, 0; t61, t62, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t65 = qJ(1) + pkin(9);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t64, -t63, 0, cos(qJ(1)) * pkin(1) + 0; t63, t64, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t70 = cos(qJ(3));
	t69 = sin(qJ(3));
	t68 = qJ(1) + pkin(9);
	t67 = cos(t68);
	t66 = sin(t68);
	t1 = [t67 * t70, -t67 * t69, t66, t67 * pkin(2) + t66 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t66 * t70, -t66 * t69, -t67, t66 * pkin(2) - t67 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t69, t70, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t78 = -qJ(4) - pkin(7);
	t77 = qJ(1) + pkin(9);
	t76 = qJ(3) + pkin(10);
	t75 = cos(t77);
	t74 = cos(t76);
	t73 = sin(t77);
	t72 = sin(t76);
	t71 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t75 * t74, -t75 * t72, t73, t75 * t71 - t73 * t78 + cos(qJ(1)) * pkin(1) + 0; t73 * t74, -t73 * t72, -t75, t73 * t71 + t75 * t78 + sin(qJ(1)) * pkin(1) + 0; t72, t74, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->23), mult. (35->24), div. (0->0), fcn. (48->10), ass. (0->15)
	t85 = qJ(1) + pkin(9);
	t81 = sin(t85);
	t87 = sin(qJ(5));
	t93 = t81 * t87;
	t88 = cos(qJ(5));
	t92 = t81 * t88;
	t83 = cos(t85);
	t91 = t83 * t87;
	t90 = t83 * t88;
	t84 = qJ(3) + pkin(10);
	t80 = sin(t84);
	t82 = cos(t84);
	t89 = pkin(4) * t82 + pkin(8) * t80 + cos(qJ(3)) * pkin(3) + pkin(2);
	t86 = -qJ(4) - pkin(7);
	t1 = [t82 * t90 + t93, -t82 * t91 + t92, t83 * t80, cos(qJ(1)) * pkin(1) - t81 * t86 + 0 + t89 * t83; t82 * t92 - t91, -t82 * t93 - t90, t81 * t80, sin(qJ(1)) * pkin(1) + t83 * t86 + 0 + t89 * t81; t80 * t88, -t80 * t87, -t82, t80 * pkin(4) - t82 * pkin(8) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:08
	% EndTime: 2020-11-04 21:36:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->28), mult. (55->30), div. (0->0), fcn. (72->10), ass. (0->19)
	t104 = qJ(1) + pkin(9);
	t100 = sin(t104);
	t106 = sin(qJ(5));
	t112 = t100 * t106;
	t107 = cos(qJ(5));
	t111 = t100 * t107;
	t102 = cos(t104);
	t110 = t102 * t106;
	t109 = t102 * t107;
	t103 = qJ(3) + pkin(10);
	t101 = cos(t103);
	t99 = sin(t103);
	t108 = pkin(4) * t101 + pkin(8) * t99 + cos(qJ(3)) * pkin(3) + pkin(2);
	t105 = -qJ(4) - pkin(7);
	t97 = t101 * t109 + t112;
	t96 = t101 * t110 - t111;
	t95 = t101 * t111 - t110;
	t94 = t101 * t112 + t109;
	t1 = [t97, t102 * t99, t96, cos(qJ(1)) * pkin(1) + t97 * pkin(5) + t96 * qJ(6) - t100 * t105 + 0 + t108 * t102; t95, t100 * t99, t94, sin(qJ(1)) * pkin(1) + t95 * pkin(5) + t94 * qJ(6) + t102 * t105 + 0 + t108 * t100; t99 * t107, -t101, t99 * t106, sin(qJ(3)) * pkin(3) - t101 * pkin(8) + pkin(6) + qJ(2) + 0 + (pkin(5) * t107 + qJ(6) * t106 + pkin(4)) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
end