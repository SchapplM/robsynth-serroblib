% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:40
	% EndTime: 2020-11-04 20:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:40
	% EndTime: 2020-11-04 20:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:40
	% EndTime: 2020-11-04 20:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t51 = cos(qJ(1));
	t50 = cos(qJ(2));
	t49 = sin(qJ(1));
	t48 = sin(qJ(2));
	t1 = [t51 * t50, -t51 * t48, t49, t51 * pkin(1) + t49 * pkin(6) + 0; t49 * t50, -t49 * t48, -t51, t49 * pkin(1) - t51 * pkin(6) + 0; t48, t50, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:40
	% EndTime: 2020-11-04 20:46:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t55 = sin(qJ(1));
	t57 = cos(qJ(2));
	t61 = t55 * t57;
	t53 = sin(qJ(3));
	t58 = cos(qJ(1));
	t60 = t58 * t53;
	t56 = cos(qJ(3));
	t59 = t58 * t56;
	t54 = sin(qJ(2));
	t52 = t57 * pkin(2) + t54 * pkin(7) + pkin(1);
	t1 = [t55 * t53 + t57 * t59, t55 * t56 - t57 * t60, t58 * t54, t55 * pkin(6) + t52 * t58 + 0; t56 * t61 - t60, -t53 * t61 - t59, t55 * t54, -t58 * pkin(6) + t52 * t55 + 0; t54 * t56, -t54 * t53, -t57, t54 * pkin(2) - t57 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:41
	% EndTime: 2020-11-04 20:46:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t69 = sin(qJ(1));
	t70 = cos(qJ(2));
	t75 = t69 * t70;
	t67 = qJ(3) + qJ(4);
	t65 = sin(t67);
	t71 = cos(qJ(1));
	t74 = t71 * t65;
	t66 = cos(t67);
	t73 = t71 * t66;
	t72 = pkin(8) + pkin(7);
	t68 = sin(qJ(2));
	t64 = cos(qJ(3)) * pkin(3) + pkin(2);
	t63 = sin(qJ(3)) * pkin(3) + pkin(6);
	t62 = t64 * t70 + t72 * t68 + pkin(1);
	t1 = [t69 * t65 + t70 * t73, t69 * t66 - t70 * t74, t71 * t68, t62 * t71 + t63 * t69 + 0; t66 * t75 - t74, -t65 * t75 - t73, t69 * t68, t62 * t69 - t63 * t71 + 0; t68 * t66, -t68 * t65, -t70, t68 * t64 - t70 * t72 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:41
	% EndTime: 2020-11-04 20:46:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->21), mult. (42->25), div. (0->0), fcn. (55->8), ass. (0->15)
	t81 = qJ(3) + qJ(4);
	t78 = sin(t81);
	t90 = pkin(6) + pkin(4) * t78 + sin(qJ(3)) * pkin(3);
	t83 = sin(qJ(1));
	t84 = cos(qJ(2));
	t89 = t83 * t84;
	t85 = cos(qJ(1));
	t88 = t85 * t78;
	t79 = cos(t81);
	t87 = t85 * t79;
	t76 = pkin(4) * t79 + cos(qJ(3)) * pkin(3) + pkin(2);
	t80 = -qJ(5) - pkin(8) - pkin(7);
	t82 = sin(qJ(2));
	t86 = t76 * t84 - t80 * t82 + pkin(1);
	t1 = [t83 * t78 + t84 * t87, t83 * t79 - t84 * t88, t85 * t82, t90 * t83 + t86 * t85 + 0; t79 * t89 - t88, -t78 * t89 - t87, t83 * t82, t86 * t83 - t90 * t85 + 0; t82 * t79, -t82 * t78, -t84, t82 * t76 + t84 * t80 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end