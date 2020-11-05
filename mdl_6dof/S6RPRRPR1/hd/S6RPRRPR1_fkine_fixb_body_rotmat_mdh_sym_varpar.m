% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t59 = qJ(1) + pkin(10);
	t58 = cos(t59);
	t57 = sin(t59);
	t1 = [t58, -t57, 0, cos(qJ(1)) * pkin(1) + 0; t57, t58, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t64 = cos(qJ(3));
	t63 = sin(qJ(3));
	t62 = qJ(1) + pkin(10);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t61 * t64, -t61 * t63, t60, t61 * pkin(2) + t60 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t60 * t64, -t60 * t63, -t61, t60 * pkin(2) - t61 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t63, t64, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t72 = -pkin(8) - pkin(7);
	t71 = qJ(3) + qJ(4);
	t70 = qJ(1) + pkin(10);
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
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->19), mult. (16->14), div. (0->0), fcn. (24->10), ass. (0->10)
	t81 = qJ(3) + qJ(4);
	t80 = qJ(1) + pkin(10);
	t79 = -qJ(5) - pkin(8) - pkin(7);
	t78 = pkin(11) + t81;
	t77 = cos(t80);
	t76 = sin(t80);
	t75 = cos(t78);
	t74 = sin(t78);
	t73 = pkin(4) * cos(t81) + cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t77 * t75, -t77 * t74, t76, t77 * t73 - t76 * t79 + cos(qJ(1)) * pkin(1) + 0; t76 * t75, -t76 * t74, -t77, t76 * t73 + t77 * t79 + sin(qJ(1)) * pkin(1) + 0; t74, t75, 0, pkin(4) * sin(t81) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:46:35
	% EndTime: 2020-11-04 21:46:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->26), mult. (38->26), div. (0->0), fcn. (51->12), ass. (0->16)
	t89 = qJ(1) + pkin(10);
	t85 = sin(t89);
	t91 = sin(qJ(6));
	t97 = t85 * t91;
	t92 = cos(qJ(6));
	t96 = t85 * t92;
	t86 = cos(t89);
	t95 = t86 * t91;
	t94 = t86 * t92;
	t90 = qJ(3) + qJ(4);
	t87 = pkin(11) + t90;
	t83 = sin(t87);
	t84 = cos(t87);
	t93 = pkin(5) * t84 + pkin(9) * t83 + pkin(4) * cos(t90) + cos(qJ(3)) * pkin(3) + pkin(2);
	t88 = -qJ(5) - pkin(8) - pkin(7);
	t1 = [t84 * t94 + t97, -t84 * t95 + t96, t86 * t83, cos(qJ(1)) * pkin(1) - t85 * t88 + 0 + t93 * t86; t84 * t96 - t95, -t84 * t97 - t94, t85 * t83, sin(qJ(1)) * pkin(1) + t86 * t88 + 0 + t93 * t85; t83 * t92, -t83 * t91, -t84, t83 * pkin(5) - t84 * pkin(9) + pkin(4) * sin(t90) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end