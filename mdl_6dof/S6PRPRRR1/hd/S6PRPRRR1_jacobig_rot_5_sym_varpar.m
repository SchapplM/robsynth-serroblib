% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->5), mult. (42->12), div. (0->0), fcn. (65->8), ass. (0->15)
t90 = sin(pkin(12));
t93 = cos(pkin(12));
t96 = sin(qJ(2));
t97 = cos(qJ(2));
t98 = t90 * t96 - t93 * t97;
t95 = cos(pkin(6));
t94 = cos(pkin(11));
t92 = sin(pkin(6));
t91 = sin(pkin(11));
t89 = -t97 * t90 - t96 * t93;
t88 = t98 * t95;
t87 = t98 * t92;
t86 = -t91 * t88 - t94 * t89;
t85 = t94 * t88 - t91 * t89;
t1 = [0, t91 * t92, 0, t86, t86, 0; 0, -t94 * t92, 0, t85, t85, 0; 0, t95, 0, t87, t87, 0;];
Jg_rot  = t1;