% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t59 = qJ(1) + qJ(2);
	t58 = cos(t59);
	t57 = sin(t59);
	t1 = [t58, -t57, 0, cos(qJ(1)) * pkin(1) + 0; t57, t58, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t63 = qJ(1) + qJ(2);
	t62 = qJ(3) + t63;
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [t61, -t60, 0, pkin(2) * cos(t63) + cos(qJ(1)) * pkin(1) + 0; t60, t61, 0, pkin(2) * sin(t63) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(8) + pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (47->19), mult. (27->23), div. (0->0), fcn. (40->10), ass. (0->13)
	t67 = qJ(1) + qJ(2);
	t66 = qJ(3) + t67;
	t64 = sin(t66);
	t68 = sin(pkin(5));
	t75 = t64 * t68;
	t65 = cos(t66);
	t74 = t65 * t68;
	t69 = cos(pkin(5));
	t70 = sin(qJ(4));
	t73 = t69 * t70;
	t71 = cos(qJ(4));
	t72 = t69 * t71;
	t1 = [-t64 * t73 + t65 * t71, -t64 * t72 - t65 * t70, t75, t65 * pkin(3) + pkin(9) * t75 + pkin(2) * cos(t67) + cos(qJ(1)) * pkin(1) + 0; t64 * t71 + t65 * t73, -t64 * t70 + t65 * t72, -t74, t64 * pkin(3) - pkin(9) * t74 + pkin(2) * sin(t67) + sin(qJ(1)) * pkin(1) + 0; t68 * t70, t68 * t71, t69, t69 * pkin(9) + pkin(6) + pkin(7) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:30:11
	% EndTime: 2024-09-27 17:30:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (89->30), mult. (41->30), div. (0->0), fcn. (48->16), ass. (0->22)
	t97 = pkin(4) * sin(qJ(4));
	t92 = qJ(1) + qJ(2);
	t91 = qJ(4) + qJ(5);
	t96 = pkin(9) + pkin(10);
	t94 = cos(pkin(5));
	t93 = sin(pkin(5));
	t90 = qJ(3) + t92;
	t89 = cos(t91);
	t88 = sin(t91);
	t87 = pkin(5) - t91;
	t86 = pkin(5) + t91;
	t85 = cos(qJ(4)) * pkin(4) + pkin(3);
	t84 = cos(t90);
	t83 = sin(t90);
	t82 = cos(t86);
	t81 = sin(t87);
	t80 = cos(t87) / 0.2e1;
	t79 = sin(t86) / 0.2e1;
	t78 = -t93 * t96 + t94 * t97;
	t77 = t80 + t82 / 0.2e1;
	t76 = t79 - t81 / 0.2e1;
	t1 = [-t83 * t76 + t84 * t89, -t83 * t77 - t84 * t88, t83 * t93, t84 * t85 - t83 * t78 + pkin(2) * cos(t92) + cos(qJ(1)) * pkin(1) + 0; t84 * t76 + t83 * t89, t84 * t77 - t83 * t88, -t84 * t93, t83 * t85 + t84 * t78 + pkin(2) * sin(t92) + sin(qJ(1)) * pkin(1) + 0; t80 - t82 / 0.2e1, t79 + t81 / 0.2e1, t94, t93 * t97 + t94 * t96 + pkin(6) + pkin(7) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end