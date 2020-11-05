% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP6 (for one body)
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

function Tc_mdh = S5PRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t58 = cos(qJ(2));
	t57 = sin(qJ(2));
	t56 = cos(pkin(8));
	t55 = sin(pkin(8));
	t1 = [t56 * t58, -t56 * t57, t55, t56 * pkin(1) + t55 * pkin(5) + 0; t55 * t58, -t55 * t57, -t56, t55 * pkin(1) - t56 * pkin(5) + 0; t57, t58, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t61 = sin(qJ(3));
	t64 = cos(qJ(2));
	t67 = t61 * t64;
	t63 = cos(qJ(3));
	t66 = t63 * t64;
	t62 = sin(qJ(2));
	t65 = pkin(2) * t64 + pkin(6) * t62 + pkin(1);
	t60 = cos(pkin(8));
	t59 = sin(pkin(8));
	t1 = [t59 * t61 + t60 * t66, t59 * t63 - t60 * t67, t60 * t62, t59 * pkin(5) + t65 * t60 + 0; t59 * t66 - t60 * t61, -t59 * t67 - t60 * t63, t59 * t62, -t60 * pkin(5) + t65 * t59 + 0; t62 * t63, -t62 * t61, -t64, t62 * pkin(2) - t64 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (37->24), div. (0->0), fcn. (50->8), ass. (0->14)
	t72 = sin(pkin(8));
	t76 = cos(qJ(2));
	t81 = t72 * t76;
	t73 = cos(pkin(8));
	t80 = t73 * t76;
	t79 = pkin(3) * sin(qJ(3)) + pkin(5);
	t68 = cos(qJ(3)) * pkin(3) + pkin(2);
	t75 = sin(qJ(2));
	t77 = pkin(7) + pkin(6);
	t78 = t68 * t76 + t75 * t77 + pkin(1);
	t71 = qJ(3) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t1 = [t72 * t69 + t70 * t80, -t69 * t80 + t72 * t70, t73 * t75, t79 * t72 + t78 * t73 + 0; -t73 * t69 + t70 * t81, -t69 * t81 - t73 * t70, t72 * t75, t78 * t72 - t79 * t73 + 0; t75 * t70, -t75 * t69, -t76, t75 * t68 - t76 * t77 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:41
	% EndTime: 2020-11-04 20:06:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->24), mult. (57->30), div. (0->0), fcn. (74->8), ass. (0->18)
	t90 = sin(pkin(8));
	t94 = cos(qJ(2));
	t99 = t90 * t94;
	t91 = cos(pkin(8));
	t98 = t91 * t94;
	t97 = pkin(3) * sin(qJ(3)) + pkin(5);
	t86 = cos(qJ(3)) * pkin(3) + pkin(2);
	t93 = sin(qJ(2));
	t95 = pkin(7) + pkin(6);
	t96 = t86 * t94 + t93 * t95 + pkin(1);
	t89 = qJ(3) + qJ(4);
	t88 = cos(t89);
	t87 = sin(t89);
	t85 = t90 * t87 + t88 * t98;
	t84 = t87 * t98 - t90 * t88;
	t83 = -t91 * t87 + t88 * t99;
	t82 = t87 * t99 + t91 * t88;
	t1 = [t85, t91 * t93, t84, t85 * pkin(4) + t84 * qJ(5) + t97 * t90 + t96 * t91 + 0; t83, t90 * t93, t82, t83 * pkin(4) + t82 * qJ(5) + t96 * t90 - t97 * t91 + 0; t93 * t88, -t94, t93 * t87, -t94 * t95 + qJ(1) + 0 + (pkin(4) * t88 + qJ(5) * t87 + t86) * t93; 0, 0, 0, 1;];
	Tc_mdh = t1;
end