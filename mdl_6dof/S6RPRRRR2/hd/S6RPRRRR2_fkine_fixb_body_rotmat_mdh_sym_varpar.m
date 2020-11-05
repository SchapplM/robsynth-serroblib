% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR2 (for one body)
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
% Datum: 2020-11-04 21:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t59 = qJ(1) + pkin(11);
	t58 = cos(t59);
	t57 = sin(t59);
	t1 = [t58, -t57, 0, cos(qJ(1)) * pkin(1) + 0; t57, t58, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t64 = cos(qJ(3));
	t63 = sin(qJ(3));
	t62 = qJ(1) + pkin(11);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t61 * t64, -t61 * t63, t60, t61 * pkin(2) + t60 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t60 * t64, -t60 * t63, -t61, t60 * pkin(2) - t61 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t63, t64, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t72 = -pkin(8) - pkin(7);
	t71 = qJ(3) + qJ(4);
	t70 = qJ(1) + pkin(11);
	t69 = cos(t71);
	t68 = sin(t71);
	t67 = cos(t70);
	t66 = sin(t70);
	t65 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t67 * t69, -t67 * t68, t66, t67 * t65 - t66 * t72 + cos(qJ(1)) * pkin(1) + 0; t66 * t69, -t66 * t68, -t67, t66 * t65 + t67 * t72 + sin(qJ(1)) * pkin(1) + 0; t68, t69, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->23), mult. (35->26), div. (0->0), fcn. (48->10), ass. (0->13)
	t79 = qJ(3) + qJ(4);
	t77 = cos(t79);
	t80 = sin(qJ(5));
	t85 = t77 * t80;
	t81 = cos(qJ(5));
	t84 = t77 * t81;
	t76 = sin(t79);
	t83 = pkin(4) * t77 + pkin(9) * t76 + cos(qJ(3)) * pkin(3) + pkin(2);
	t82 = -pkin(8) - pkin(7);
	t78 = qJ(1) + pkin(11);
	t75 = cos(t78);
	t74 = sin(t78);
	t1 = [t74 * t80 + t75 * t84, t74 * t81 - t75 * t85, t75 * t76, cos(qJ(1)) * pkin(1) - t74 * t82 + 0 + t83 * t75; t74 * t84 - t75 * t80, -t74 * t85 - t75 * t81, t74 * t76, sin(qJ(1)) * pkin(1) + t75 * t82 + 0 + t83 * t74; t76 * t81, -t76 * t80, -t77, t76 * pkin(4) - t77 * pkin(9) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:55:05
	% EndTime: 2020-11-04 21:55:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->27), mult. (42->28), div. (0->0), fcn. (55->12), ass. (0->16)
	t95 = qJ(5) + qJ(6);
	t90 = sin(t95);
	t96 = qJ(3) + qJ(4);
	t93 = cos(t96);
	t103 = t90 * t93;
	t92 = cos(t95);
	t102 = t92 * t93;
	t101 = pkin(5) * sin(qJ(5)) + pkin(8) + pkin(7);
	t86 = cos(qJ(5)) * pkin(5) + pkin(4);
	t91 = sin(t96);
	t98 = -pkin(10) - pkin(9);
	t100 = t86 * t93 - t91 * t98 + cos(qJ(3)) * pkin(3) + pkin(2);
	t94 = qJ(1) + pkin(11);
	t89 = cos(t94);
	t88 = sin(t94);
	t1 = [t89 * t102 + t88 * t90, -t89 * t103 + t88 * t92, t89 * t91, cos(qJ(1)) * pkin(1) + 0 + t101 * t88 + t100 * t89; t88 * t102 - t89 * t90, -t88 * t103 - t89 * t92, t88 * t91, sin(qJ(1)) * pkin(1) + 0 - t101 * t89 + t100 * t88; t91 * t92, -t91 * t90, -t93, t91 * t86 + t93 * t98 + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end