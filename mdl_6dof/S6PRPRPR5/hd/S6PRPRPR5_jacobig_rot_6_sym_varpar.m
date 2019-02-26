% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t94 = sin(pkin(10));
t95 = sin(pkin(6));
t103 = t94 * t95;
t96 = cos(pkin(10));
t102 = t96 * t95;
t97 = cos(pkin(6));
t98 = sin(qJ(2));
t101 = t97 * t98;
t99 = cos(qJ(2));
t100 = t97 * t99;
t93 = pkin(11) + qJ(4);
t92 = cos(t93);
t91 = sin(t93);
t1 = [0, t103, 0, t94 * t100 + t96 * t98, 0 (-t94 * t101 + t96 * t99) * t92 + t91 * t103; 0, -t102, 0, -t96 * t100 + t94 * t98, 0 (t96 * t101 + t94 * t99) * t92 - t91 * t102; 0, t97, 0, -t95 * t99, 0, t95 * t98 * t92 + t97 * t91;];
Jg_rot  = t1;
