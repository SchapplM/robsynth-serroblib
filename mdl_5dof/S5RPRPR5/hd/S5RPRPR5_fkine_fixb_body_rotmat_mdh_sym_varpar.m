% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t56 = cos(pkin(8));
	t55 = sin(pkin(8));
	t1 = [t58 * t56, -t58 * t55, t57, t58 * pkin(1) + t57 * qJ(2) + 0; t57 * t56, -t57 * t55, -t58, t57 * pkin(1) - t58 * qJ(2) + 0; t55, t56, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t62 = sin(qJ(3));
	t63 = sin(qJ(1));
	t69 = t63 * t62;
	t64 = cos(qJ(3));
	t68 = t63 * t64;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t61 = cos(pkin(8));
	t60 = sin(pkin(8));
	t59 = pkin(2) * t61 + t60 * pkin(6) + pkin(1);
	t1 = [t61 * t66 + t69, -t61 * t67 + t68, t65 * t60, t63 * qJ(2) + t59 * t65 + 0; t61 * t68 - t67, -t61 * t69 - t66, t63 * t60, -t65 * qJ(2) + t59 * t63 + 0; t60 * t64, -t60 * t62, -t61, t60 * pkin(2) - t61 * pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (34->22), div. (0->0), fcn. (47->8), ass. (0->16)
	t74 = qJ(3) + pkin(9);
	t72 = sin(t74);
	t78 = sin(qJ(1));
	t85 = t78 * t72;
	t73 = cos(t74);
	t84 = t78 * t73;
	t80 = cos(qJ(1));
	t83 = t80 * t72;
	t82 = t80 * t73;
	t81 = pkin(3) * cos(qJ(3)) + pkin(2);
	t77 = qJ(4) + pkin(6);
	t76 = cos(pkin(8));
	t75 = sin(pkin(8));
	t71 = sin(qJ(3)) * pkin(3) + qJ(2);
	t70 = t77 * t75 + t81 * t76 + pkin(1);
	t1 = [t76 * t82 + t85, -t76 * t83 + t84, t80 * t75, t70 * t80 + t71 * t78 + 0; t76 * t84 - t83, -t76 * t85 - t82, t78 * t75, t70 * t78 - t71 * t80 + 0; t75 * t73, -t75 * t72, -t76, t81 * t75 - t76 * t77 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:24:45
	% EndTime: 2022-01-23 09:24:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t92 = qJ(3) + pkin(9);
	t90 = qJ(5) + t92;
	t88 = sin(t90);
	t95 = sin(qJ(1));
	t102 = t95 * t88;
	t89 = cos(t90);
	t101 = t95 * t89;
	t96 = cos(qJ(1));
	t100 = t96 * t88;
	t99 = t96 * t89;
	t98 = qJ(2) + pkin(4) * sin(t92) + sin(qJ(3)) * pkin(3);
	t86 = pkin(4) * cos(t92) + pkin(3) * cos(qJ(3)) + pkin(2);
	t91 = -pkin(7) - qJ(4) - pkin(6);
	t93 = sin(pkin(8));
	t94 = cos(pkin(8));
	t97 = t86 * t94 - t91 * t93 + pkin(1);
	t1 = [t94 * t99 + t102, -t94 * t100 + t101, t96 * t93, t98 * t95 + t97 * t96 + 0; t94 * t101 - t100, -t94 * t102 - t99, t95 * t93, t97 * t95 - t98 * t96 + 0; t93 * t89, -t93 * t88, -t94, t93 * t86 + t94 * t91 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end