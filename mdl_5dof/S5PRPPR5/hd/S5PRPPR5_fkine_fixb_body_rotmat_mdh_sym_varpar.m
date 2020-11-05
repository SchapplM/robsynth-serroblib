% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:57
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t55 = cos(pkin(7));
	t54 = sin(pkin(7));
	t1 = [t55, -t54, 0, 0; t54, t55, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t59 = cos(qJ(2));
	t58 = sin(qJ(2));
	t57 = cos(pkin(7));
	t56 = sin(pkin(7));
	t1 = [t57 * t59, -t57 * t58, t56, t57 * pkin(1) + t56 * pkin(5) + 0; t56 * t59, -t56 * t58, -t57, t56 * pkin(1) - t57 * pkin(5) + 0; t58, t59, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (18->12), div. (0->0), fcn. (26->4), ass. (0->6)
	t62 = sin(qJ(2));
	t63 = cos(qJ(2));
	t64 = pkin(2) * t63 + qJ(3) * t62 + pkin(1);
	t61 = cos(pkin(7));
	t60 = sin(pkin(7));
	t1 = [t61 * t63, t60, t61 * t62, t60 * pkin(5) + t64 * t61 + 0; t60 * t63, -t61, t60 * t62, -t61 * pkin(5) + t64 * t60 + 0; t62, 0, -t63, t62 * pkin(2) - t63 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (30->16), div. (0->0), fcn. (44->6), ass. (0->12)
	t77 = cos(pkin(8));
	t76 = sin(pkin(8));
	t72 = sin(qJ(2));
	t73 = cos(qJ(2));
	t74 = pkin(2) + pkin(3);
	t75 = qJ(3) * t72 + t73 * t74 + pkin(1);
	t71 = pkin(5) - qJ(4);
	t70 = cos(pkin(7));
	t69 = sin(pkin(7));
	t66 = t72 * t77 - t73 * t76;
	t65 = -t72 * t76 - t73 * t77;
	t1 = [-t70 * t65, t70 * t66, -t69, t71 * t69 + t75 * t70 + 0; -t69 * t65, t69 * t66, t70, t75 * t69 - t71 * t70 + 0; t66, t65, 0, -t73 * qJ(3) + t74 * t72 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:57:31
	% EndTime: 2020-11-04 19:57:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->24), mult. (60->28), div. (0->0), fcn. (82->8), ass. (0->19)
	t84 = sin(pkin(7));
	t88 = sin(qJ(5));
	t96 = t84 * t88;
	t90 = cos(qJ(5));
	t95 = t84 * t90;
	t86 = cos(pkin(7));
	t94 = t86 * t88;
	t93 = t86 * t90;
	t83 = sin(pkin(8));
	t85 = cos(pkin(8));
	t80 = t85 * pkin(4) + t83 * pkin(6) + pkin(2) + pkin(3);
	t81 = -t83 * pkin(4) + t85 * pkin(6) - qJ(3);
	t89 = sin(qJ(2));
	t91 = cos(qJ(2));
	t92 = t80 * t91 - t81 * t89 + pkin(1);
	t87 = pkin(5) - qJ(4);
	t79 = t89 * t83 + t91 * t85;
	t78 = -t91 * t83 + t89 * t85;
	t1 = [t79 * t93 - t96, -t79 * t94 - t95, -t86 * t78, t87 * t84 + t92 * t86 + 0; t79 * t95 + t94, -t79 * t96 + t93, -t84 * t78, t92 * t84 - t87 * t86 + 0; t78 * t90, -t78 * t88, t79, t80 * t89 + t81 * t91 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end