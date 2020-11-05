% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t55 = qJ(1) + pkin(10);
	t54 = cos(t55);
	t53 = sin(t55);
	t1 = [t54, -t53, 0, cos(qJ(1)) * pkin(1) + 0; t53, t54, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t58 = qJ(1) + pkin(10);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [0, -t57, t56, t57 * pkin(2) + t56 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t56, -t57, t56 * pkin(2) - t57 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t64 = pkin(2) + pkin(7);
	t63 = cos(qJ(4));
	t62 = sin(qJ(4));
	t61 = qJ(1) + pkin(10);
	t60 = cos(t61);
	t59 = sin(t61);
	t1 = [t59 * t62, t59 * t63, t60, t64 * t60 + t59 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t60 * t62, -t60 * t63, t59, t64 * t59 - t60 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t63, -t62, 0, pkin(3) + qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->21), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->12)
	t68 = sin(qJ(5));
	t69 = sin(qJ(4));
	t75 = t68 * t69;
	t70 = cos(qJ(5));
	t74 = t69 * t70;
	t71 = cos(qJ(4));
	t73 = pkin(4) * t69 - pkin(8) * t71 + qJ(3);
	t72 = pkin(2) + pkin(7);
	t67 = qJ(1) + pkin(10);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [t65 * t74 + t66 * t68, -t65 * t75 + t66 * t70, -t65 * t71, cos(qJ(1)) * pkin(1) + t72 * t66 + 0 + t73 * t65; t65 * t68 - t66 * t74, t65 * t70 + t66 * t75, t66 * t71, sin(qJ(1)) * pkin(1) + t72 * t65 + 0 - t73 * t66; t71 * t70, -t71 * t68, t69, t71 * pkin(4) + t69 * pkin(8) + pkin(3) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:30:56
	% EndTime: 2020-11-04 21:30:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->25), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t82 = qJ(5) + qJ(6);
	t79 = sin(t82);
	t84 = sin(qJ(4));
	t91 = t79 * t84;
	t80 = cos(t82);
	t90 = t80 * t84;
	t89 = pkin(5) * sin(qJ(5)) + pkin(2) + pkin(7);
	t76 = cos(qJ(5)) * pkin(5) + pkin(4);
	t85 = cos(qJ(4));
	t86 = -pkin(9) - pkin(8);
	t88 = t76 * t84 + t85 * t86 + qJ(3);
	t81 = qJ(1) + pkin(10);
	t78 = cos(t81);
	t77 = sin(t81);
	t1 = [t77 * t90 + t78 * t79, -t77 * t91 + t78 * t80, -t77 * t85, cos(qJ(1)) * pkin(1) + 0 + t89 * t78 + t88 * t77; t77 * t79 - t78 * t90, t77 * t80 + t78 * t91, t78 * t85, sin(qJ(1)) * pkin(1) + 0 + t89 * t77 - t88 * t78; t85 * t80, -t85 * t79, t84, t85 * t76 - t84 * t86 + pkin(3) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end