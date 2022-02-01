% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t1 = [t55, -t54, 0, 0; t54, t55, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t57 = cos(pkin(8));
	t56 = sin(pkin(8));
	t1 = [t59 * t57, -t59 * t56, t58, t59 * pkin(1) + t58 * qJ(2) + 0; t58 * t57, -t58 * t56, -t59, t58 * pkin(1) - t59 * qJ(2) + 0; t56, t57, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t63 = sin(qJ(3));
	t64 = sin(qJ(1));
	t70 = t64 * t63;
	t65 = cos(qJ(3));
	t69 = t64 * t65;
	t66 = cos(qJ(1));
	t68 = t66 * t63;
	t67 = t66 * t65;
	t62 = cos(pkin(8));
	t61 = sin(pkin(8));
	t60 = pkin(2) * t62 + t61 * pkin(6) + pkin(1);
	t1 = [t62 * t67 + t70, -t62 * t68 + t69, t66 * t61, t64 * qJ(2) + t60 * t66 + 0; t62 * t69 - t68, -t62 * t70 - t67, t64 * t61, -t66 * qJ(2) + t60 * t64 + 0; t61 * t65, -t61 * t63, -t62, t61 * pkin(2) - t62 * pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (34->22), div. (0->0), fcn. (47->8), ass. (0->16)
	t75 = qJ(3) + qJ(4);
	t73 = sin(t75);
	t78 = sin(qJ(1));
	t86 = t78 * t73;
	t74 = cos(t75);
	t85 = t78 * t74;
	t80 = cos(qJ(1));
	t84 = t80 * t73;
	t83 = t80 * t74;
	t82 = pkin(3) * cos(qJ(3)) + pkin(2);
	t81 = pkin(7) + pkin(6);
	t77 = cos(pkin(8));
	t76 = sin(pkin(8));
	t72 = sin(qJ(3)) * pkin(3) + qJ(2);
	t71 = t81 * t76 + t82 * t77 + pkin(1);
	t1 = [t77 * t83 + t86, -t77 * t84 + t85, t80 * t76, t71 * t80 + t72 * t78 + 0; t77 * t85 - t84, -t77 * t86 - t83, t78 * t76, t71 * t78 - t72 * t80 + 0; t76 * t74, -t76 * t73, -t77, t82 * t76 - t77 * t81 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:31:40
	% EndTime: 2022-01-23 09:31:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->21), mult. (42->24), div. (0->0), fcn. (55->8), ass. (0->16)
	t92 = qJ(3) + qJ(4);
	t89 = sin(t92);
	t95 = sin(qJ(1));
	t102 = t95 * t89;
	t90 = cos(t92);
	t101 = t95 * t90;
	t96 = cos(qJ(1));
	t100 = t96 * t89;
	t99 = t96 * t90;
	t98 = qJ(2) + pkin(4) * t89 + sin(qJ(3)) * pkin(3);
	t87 = pkin(4) * t90 + pkin(3) * cos(qJ(3)) + pkin(2);
	t91 = -qJ(5) - pkin(7) - pkin(6);
	t93 = sin(pkin(8));
	t94 = cos(pkin(8));
	t97 = t87 * t94 - t91 * t93 + pkin(1);
	t1 = [t94 * t99 + t102, -t94 * t100 + t101, t96 * t93, t98 * t95 + t97 * t96 + 0; t94 * t101 - t100, -t94 * t102 - t99, t95 * t93, t97 * t95 - t98 * t96 + 0; t93 * t90, -t93 * t89, -t94, t93 * t87 + t94 * t91 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end