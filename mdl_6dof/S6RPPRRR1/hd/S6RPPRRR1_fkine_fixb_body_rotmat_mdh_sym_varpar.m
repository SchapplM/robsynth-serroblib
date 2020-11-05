% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t60 = qJ(1) + pkin(10);
	t59 = cos(t60);
	t58 = sin(t60);
	t1 = [t59, -t58, 0, cos(qJ(1)) * pkin(1) + 0; t58, t59, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t65 = cos(pkin(11));
	t64 = sin(pkin(11));
	t63 = qJ(1) + pkin(10);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62 * t65, -t62 * t64, t61, t62 * pkin(2) + t61 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; t61 * t65, -t61 * t64, -t62, t61 * pkin(2) - t62 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t64, t65, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t73 = -pkin(7) - qJ(3);
	t72 = qJ(1) + pkin(10);
	t71 = pkin(11) + qJ(4);
	t70 = cos(t72);
	t69 = cos(t71);
	t68 = sin(t72);
	t67 = sin(t71);
	t66 = cos(pkin(11)) * pkin(3) + pkin(2);
	t1 = [t70 * t69, -t70 * t67, t68, t70 * t66 - t68 * t73 + cos(qJ(1)) * pkin(1) + 0; t68 * t69, -t68 * t67, -t70, t68 * t66 + t70 * t73 + sin(qJ(1)) * pkin(1) + 0; t67, t69, 0, sin(pkin(11)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->19), mult. (16->14), div. (0->0), fcn. (24->10), ass. (0->10)
	t81 = pkin(11) + qJ(4);
	t82 = qJ(1) + pkin(10);
	t80 = qJ(3) + pkin(7) + pkin(8);
	t79 = qJ(5) + t81;
	t78 = cos(t82);
	t77 = sin(t82);
	t76 = cos(t79);
	t75 = sin(t79);
	t74 = pkin(4) * cos(t81) + cos(pkin(11)) * pkin(3) + pkin(2);
	t1 = [t78 * t76, -t78 * t75, t77, t78 * t74 + t80 * t77 + cos(qJ(1)) * pkin(1) + 0; t77 * t76, -t77 * t75, -t78, t77 * t74 - t78 * t80 + sin(qJ(1)) * pkin(1) + 0; t75, t76, 0, pkin(4) * sin(t81) + sin(pkin(11)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:15
	% EndTime: 2020-11-04 21:30:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (81->26), mult. (38->26), div. (0->0), fcn. (51->12), ass. (0->16)
	t91 = qJ(1) + pkin(10);
	t86 = sin(t91);
	t92 = sin(qJ(6));
	t98 = t86 * t92;
	t93 = cos(qJ(6));
	t97 = t86 * t93;
	t87 = cos(t91);
	t96 = t87 * t92;
	t95 = t87 * t93;
	t90 = pkin(11) + qJ(4);
	t88 = qJ(5) + t90;
	t84 = sin(t88);
	t85 = cos(t88);
	t94 = pkin(5) * t85 + pkin(9) * t84 + pkin(4) * cos(t90) + cos(pkin(11)) * pkin(3) + pkin(2);
	t89 = qJ(3) + pkin(7) + pkin(8);
	t1 = [t85 * t95 + t98, -t85 * t96 + t97, t87 * t84, cos(qJ(1)) * pkin(1) + t89 * t86 + 0 + t94 * t87; t85 * t97 - t96, -t85 * t98 - t95, t86 * t84, sin(qJ(1)) * pkin(1) - t87 * t89 + 0 + t94 * t86; t84 * t93, -t84 * t92, -t85, t84 * pkin(5) - t85 * pkin(9) + pkin(4) * sin(t90) + sin(pkin(11)) * pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end