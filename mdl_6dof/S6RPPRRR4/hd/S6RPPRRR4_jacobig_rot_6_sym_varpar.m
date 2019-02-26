% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPPRRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->4), mult. (16->6), div. (0->0), fcn. (30->6), ass. (0->11)
t99 = cos(qJ(1));
t98 = cos(qJ(4));
t97 = sin(qJ(1));
t96 = sin(qJ(4));
t95 = cos(pkin(10));
t94 = sin(pkin(10));
t93 = t97 * t94 + t99 * t95;
t92 = t99 * t94 - t97 * t95;
t91 = t93 * t96;
t90 = t92 * t96;
t1 = [0, 0, 0, t92, t91, t91; 0, 0, 0, t93, -t90, -t90; 1, 0, 0, 0, t98, t98;];
Jg_rot  = t1;
