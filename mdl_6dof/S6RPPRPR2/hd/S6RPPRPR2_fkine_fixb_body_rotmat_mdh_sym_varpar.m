% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t60 = qJ(1) + pkin(9);
	t59 = cos(t60);
	t58 = sin(t60);
	t1 = [t59, -t58, 0, cos(qJ(1)) * pkin(1) + 0; t58, t59, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t65 = cos(pkin(10));
	t64 = sin(pkin(10));
	t63 = qJ(1) + pkin(9);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62 * t65, -t62 * t64, t61, t62 * pkin(2) + t61 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; t61 * t65, -t61 * t64, -t62, t61 * pkin(2) - t62 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t64, t65, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t73 = -pkin(7) - qJ(3);
	t72 = qJ(1) + pkin(9);
	t71 = pkin(10) + qJ(4);
	t70 = cos(t72);
	t69 = cos(t71);
	t68 = sin(t72);
	t67 = sin(t71);
	t66 = cos(pkin(10)) * pkin(3) + pkin(2);
	t1 = [t70 * t69, -t70 * t67, t68, t70 * t66 - t68 * t73 + cos(qJ(1)) * pkin(1) + 0; t68 * t69, -t68 * t67, -t70, t68 * t66 + t70 * t73 + sin(qJ(1)) * pkin(1) + 0; t67, t69, 0, sin(pkin(10)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (50->22), mult. (23->16), div. (0->0), fcn. (31->8), ass. (0->9)
	t79 = pkin(10) + qJ(4);
	t75 = sin(t79);
	t77 = cos(t79);
	t82 = pkin(4) * t77 + qJ(5) * t75 + cos(pkin(10)) * pkin(3) + pkin(2);
	t81 = -pkin(7) - qJ(3);
	t80 = qJ(1) + pkin(9);
	t78 = cos(t80);
	t76 = sin(t80);
	t1 = [t76, -t78 * t77, t78 * t75, cos(qJ(1)) * pkin(1) - t76 * t81 + 0 + t82 * t78; -t78, -t76 * t77, t76 * t75, sin(qJ(1)) * pkin(1) + t78 * t81 + 0 + t82 * t76; 0, -t75, -t77, t75 * pkin(4) - t77 * qJ(5) + sin(pkin(10)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:01
	% EndTime: 2020-11-04 21:25:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (70->24), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->16)
	t99 = pkin(4) + pkin(8);
	t98 = pkin(5) + pkin(7) + qJ(3);
	t89 = qJ(1) + pkin(9);
	t85 = sin(t89);
	t91 = sin(qJ(6));
	t97 = t85 * t91;
	t92 = cos(qJ(6));
	t96 = t85 * t92;
	t87 = cos(t89);
	t95 = t87 * t91;
	t94 = t87 * t92;
	t88 = pkin(10) + qJ(4);
	t84 = sin(t88);
	t86 = cos(t88);
	t93 = qJ(5) * t84 + t99 * t86 + cos(pkin(10)) * pkin(3) + pkin(2);
	t1 = [t84 * t95 + t96, t84 * t94 - t97, t87 * t86, cos(qJ(1)) * pkin(1) + 0 + t98 * t85 + t93 * t87; t84 * t97 - t94, t84 * t96 + t95, t85 * t86, sin(qJ(1)) * pkin(1) + 0 - t98 * t87 + t93 * t85; -t86 * t91, -t86 * t92, t84, sin(pkin(10)) * pkin(3) - t86 * qJ(5) + pkin(6) + qJ(2) + 0 + t99 * t84; 0, 0, 0, 1;];
	Tc_mdh = t1;
end