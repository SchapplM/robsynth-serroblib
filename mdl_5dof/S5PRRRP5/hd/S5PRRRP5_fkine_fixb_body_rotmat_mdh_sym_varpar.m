% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(2));
	t52 = sin(qJ(2));
	t51 = cos(pkin(8));
	t50 = sin(pkin(8));
	t1 = [t51 * t53, -t51 * t52, t50, pkin(1) * t51 + pkin(5) * t50 + 0; t50 * t53, -t50 * t52, -t51, pkin(1) * t50 - pkin(5) * t51 + 0; t52, t53, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t56 = sin(qJ(3));
	t59 = cos(qJ(2));
	t62 = t56 * t59;
	t58 = cos(qJ(3));
	t61 = t58 * t59;
	t57 = sin(qJ(2));
	t60 = pkin(2) * t59 + pkin(6) * t57 + pkin(1);
	t55 = cos(pkin(8));
	t54 = sin(pkin(8));
	t1 = [t54 * t56 + t55 * t61, t54 * t58 - t55 * t62, t55 * t57, t54 * pkin(5) + t60 * t55 + 0; t54 * t61 - t55 * t56, -t54 * t62 - t55 * t58, t54 * t57, -t55 * pkin(5) + t60 * t54 + 0; t57 * t58, -t57 * t56, -t59, t57 * pkin(2) - t59 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (35->19), mult. (37->24), div. (0->0), fcn. (50->8), ass. (0->14)
	t67 = sin(pkin(8));
	t71 = cos(qJ(2));
	t76 = t67 * t71;
	t68 = cos(pkin(8));
	t75 = t68 * t71;
	t74 = pkin(3) * sin(qJ(3)) + pkin(5);
	t63 = cos(qJ(3)) * pkin(3) + pkin(2);
	t70 = sin(qJ(2));
	t72 = pkin(7) + pkin(6);
	t73 = t63 * t71 + t70 * t72 + pkin(1);
	t66 = qJ(3) + qJ(4);
	t65 = cos(t66);
	t64 = sin(t66);
	t1 = [t67 * t64 + t65 * t75, -t64 * t75 + t67 * t65, t68 * t70, t74 * t67 + t73 * t68 + 0; -t68 * t64 + t65 * t76, -t64 * t76 - t68 * t65, t67 * t70, t73 * t67 - t74 * t68 + 0; t70 * t65, -t70 * t64, -t71, t70 * t63 - t71 * t72 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:24
	% EndTime: 2020-11-04 20:06:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->21), mult. (42->26), div. (0->0), fcn. (55->8), ass. (0->14)
	t82 = qJ(3) + qJ(4);
	t79 = sin(t82);
	t90 = pkin(5) + pkin(4) * t79 + pkin(3) * sin(qJ(3));
	t83 = sin(pkin(8));
	t86 = cos(qJ(2));
	t89 = t83 * t86;
	t84 = cos(pkin(8));
	t88 = t84 * t86;
	t80 = cos(t82);
	t77 = pkin(4) * t80 + cos(qJ(3)) * pkin(3) + pkin(2);
	t81 = -qJ(5) - pkin(7) - pkin(6);
	t85 = sin(qJ(2));
	t87 = t77 * t86 - t81 * t85 + pkin(1);
	t1 = [t83 * t79 + t80 * t88, -t79 * t88 + t83 * t80, t84 * t85, t90 * t83 + t87 * t84 + 0; -t84 * t79 + t80 * t89, -t79 * t89 - t84 * t80, t83 * t85, t87 * t83 - t90 * t84 + 0; t85 * t80, -t85 * t79, -t86, t85 * t77 + t86 * t81 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end