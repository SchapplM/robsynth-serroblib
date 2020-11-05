% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(1));
	t54 = cos(qJ(2));
	t53 = sin(qJ(1));
	t52 = sin(qJ(2));
	t1 = [t55 * t54, -t55 * t52, t53, t55 * pkin(1) + t53 * pkin(6) + 0; t53 * t54, -t53 * t52, -t55, t53 * pkin(1) - t55 * pkin(6) + 0; t52, t54, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t56 = t59 * pkin(2) + t57 * qJ(3) + pkin(1);
	t1 = [t60 * t59, t58, t60 * t57, t58 * pkin(6) + t56 * t60 + 0; t58 * t59, -t60, t58 * t57, -t60 * pkin(6) + t56 * t58 + 0; t57, 0, -t59, t57 * pkin(2) - t59 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t73 = cos(qJ(4));
	t72 = sin(qJ(4));
	t71 = pkin(2) + pkin(3);
	t70 = pkin(6) - pkin(7);
	t69 = cos(qJ(1));
	t68 = cos(qJ(2));
	t67 = sin(qJ(1));
	t66 = sin(qJ(2));
	t63 = t66 * qJ(3) + t71 * t68 + pkin(1);
	t62 = t66 * t73 - t68 * t72;
	t61 = -t66 * t72 - t68 * t73;
	t1 = [-t69 * t61, t69 * t62, -t67, t63 * t69 + t70 * t67 + 0; -t67 * t61, t67 * t62, t69, t63 * t67 - t70 * t69 + 0; t62, t61, 0, -t68 * qJ(3) + t71 * t66 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:38:14
	% EndTime: 2020-11-04 20:38:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->24), mult. (56->28), div. (0->0), fcn. (78->8), ass. (0->19)
	t80 = sin(qJ(5));
	t83 = sin(qJ(1));
	t92 = t83 * t80;
	t84 = cos(qJ(5));
	t91 = t83 * t84;
	t87 = cos(qJ(1));
	t90 = t87 * t80;
	t89 = t87 * t84;
	t88 = pkin(6) - pkin(7);
	t86 = cos(qJ(2));
	t85 = cos(qJ(4));
	t82 = sin(qJ(2));
	t81 = sin(qJ(4));
	t78 = -t81 * pkin(4) + t85 * pkin(8) - qJ(3);
	t77 = t85 * pkin(4) + t81 * pkin(8) + pkin(2) + pkin(3);
	t76 = t82 * t81 + t86 * t85;
	t75 = -t86 * t81 + t82 * t85;
	t74 = t77 * t86 - t78 * t82 + pkin(1);
	t1 = [t76 * t89 - t92, -t76 * t90 - t91, -t87 * t75, t74 * t87 + t88 * t83 + 0; t76 * t91 + t90, -t76 * t92 + t89, -t83 * t75, t74 * t83 - t88 * t87 + 0; t75 * t84, -t75 * t80, t76, t77 * t82 + t78 * t86 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end