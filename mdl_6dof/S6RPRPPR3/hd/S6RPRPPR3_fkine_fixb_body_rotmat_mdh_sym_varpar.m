% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t63 = qJ(1) + pkin(9);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62, -t61, 0, cos(qJ(1)) * pkin(1) + 0; t61, t62, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t68 = cos(qJ(3));
	t67 = sin(qJ(3));
	t66 = qJ(1) + pkin(9);
	t65 = cos(t66);
	t64 = sin(t66);
	t1 = [t65 * t68, -t65 * t67, t64, t65 * pkin(2) + t64 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t64 * t68, -t64 * t67, -t65, t64 * pkin(2) - t65 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t67, t68, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t72 = sin(qJ(3));
	t73 = cos(qJ(3));
	t74 = pkin(3) * t73 + qJ(4) * t72 + pkin(2);
	t71 = qJ(1) + pkin(9);
	t70 = cos(t71);
	t69 = sin(t71);
	t1 = [t70 * t73, t69, t70 * t72, cos(qJ(1)) * pkin(1) + t69 * pkin(7) + 0 + t74 * t70; t69 * t73, -t70, t69 * t72, sin(qJ(1)) * pkin(1) - t70 * pkin(7) + 0 + t74 * t69; t72, 0, -t73, t72 * pkin(3) - t73 * qJ(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->20), mult. (27->14), div. (0->0), fcn. (35->6), ass. (0->9)
	t82 = pkin(3) + pkin(4);
	t81 = qJ(5) - pkin(7);
	t78 = sin(qJ(3));
	t79 = cos(qJ(3));
	t80 = qJ(4) * t78 + t82 * t79 + pkin(2);
	t77 = qJ(1) + pkin(9);
	t76 = cos(t77);
	t75 = sin(t77);
	t1 = [t76 * t78, -t76 * t79, -t75, cos(qJ(1)) * pkin(1) + 0 - t81 * t75 + t80 * t76; t75 * t78, -t75 * t79, t76, sin(qJ(1)) * pkin(1) + 0 + t81 * t76 + t80 * t75; -t79, -t78, 0, -t79 * qJ(4) + t82 * t78 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:33:45
	% EndTime: 2020-11-04 21:33:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (70->29), mult. (46->32), div. (0->0), fcn. (45->12), ass. (0->19)
	t91 = qJ(1) + pkin(9);
	t89 = qJ(3) + t91;
	t102 = sin(t89) / 0.2e1;
	t90 = -qJ(3) + t91;
	t101 = cos(t90) / 0.2e1;
	t95 = sin(qJ(6));
	t96 = sin(qJ(3));
	t100 = t95 * t96;
	t97 = cos(qJ(6));
	t99 = t96 * t97;
	t98 = cos(qJ(3));
	t94 = qJ(4) + pkin(5);
	t93 = qJ(5) - pkin(7);
	t92 = pkin(3) + pkin(4) + pkin(8);
	t88 = cos(t91);
	t87 = sin(t91);
	t85 = cos(t89);
	t84 = sin(t90);
	t1 = [-t87 * t95 + t88 * t99, -t88 * t100 - t87 * t97, t88 * t98, -t93 * t87 + cos(qJ(1)) * pkin(1) + t88 * pkin(2) + 0 + (-t84 / 0.2e1 + t102) * t94 + (t101 + t85 / 0.2e1) * t92; t87 * t99 + t88 * t95, -t87 * t100 + t88 * t97, t87 * t98, t93 * t88 + sin(qJ(1)) * pkin(1) + t87 * pkin(2) + 0 + (t101 - t85 / 0.2e1) * t94 + (t84 / 0.2e1 + t102) * t92; -t98 * t97, t98 * t95, t96, t92 * t96 - t94 * t98 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end