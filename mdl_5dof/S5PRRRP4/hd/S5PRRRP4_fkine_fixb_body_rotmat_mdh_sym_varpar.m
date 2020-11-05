% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP4 (for one body)
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

function Tc_mdh = S5PRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(pkin(8));
	t51 = sin(pkin(8));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t56 = cos(qJ(2));
	t55 = sin(qJ(2));
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t1 = [t54 * t56, -t54 * t55, t53, t54 * pkin(1) + t53 * pkin(5) + 0; t53 * t56, -t53 * t55, -t54, t53 * pkin(1) - t54 * pkin(5) + 0; t55, t56, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t63 = -pkin(6) - pkin(5);
	t62 = cos(pkin(8));
	t61 = sin(pkin(8));
	t60 = qJ(2) + qJ(3);
	t59 = cos(t60);
	t58 = sin(t60);
	t57 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t62 * t59, -t62 * t58, t61, t62 * t57 - t61 * t63 + 0; t61 * t59, -t61 * t58, -t62, t61 * t57 + t62 * t63 + 0; t58, t59, 0, sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t68 = sin(pkin(8));
	t70 = sin(qJ(4));
	t77 = t68 * t70;
	t71 = cos(qJ(4));
	t76 = t68 * t71;
	t69 = cos(pkin(8));
	t75 = t69 * t70;
	t74 = t69 * t71;
	t67 = qJ(2) + qJ(3);
	t65 = sin(t67);
	t66 = cos(t67);
	t73 = pkin(3) * t66 + pkin(7) * t65 + cos(qJ(2)) * pkin(2) + pkin(1);
	t72 = -pkin(6) - pkin(5);
	t1 = [t66 * t74 + t77, -t66 * t75 + t76, t69 * t65, -t68 * t72 + t73 * t69 + 0; t66 * t76 - t75, -t66 * t77 - t74, t68 * t65, t73 * t68 + t69 * t72 + 0; t65 * t71, -t65 * t70, -t66, t65 * pkin(3) - t66 * pkin(7) + sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:08
	% EndTime: 2020-11-04 20:06:08
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t86 = sin(pkin(8));
	t88 = sin(qJ(4));
	t95 = t86 * t88;
	t89 = cos(qJ(4));
	t94 = t86 * t89;
	t87 = cos(pkin(8));
	t93 = t87 * t88;
	t92 = t87 * t89;
	t85 = qJ(2) + qJ(3);
	t83 = sin(t85);
	t84 = cos(t85);
	t91 = pkin(3) * t84 + pkin(7) * t83 + cos(qJ(2)) * pkin(2) + pkin(1);
	t90 = -pkin(6) - pkin(5);
	t81 = t84 * t92 + t95;
	t80 = t84 * t93 - t94;
	t79 = t84 * t94 - t93;
	t78 = t84 * t95 + t92;
	t1 = [t81, t87 * t83, t80, t81 * pkin(4) + t80 * qJ(5) - t86 * t90 + t91 * t87 + 0; t79, t86 * t83, t78, t79 * pkin(4) + t78 * qJ(5) + t91 * t86 + t87 * t90 + 0; t83 * t89, -t84, t83 * t88, sin(qJ(2)) * pkin(2) - t84 * pkin(7) + qJ(1) + 0 + (pkin(4) * t89 + qJ(5) * t88 + pkin(3)) * t83; 0, 0, 0, 1;];
	Tc_mdh = t1;
end