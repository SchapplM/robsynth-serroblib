% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:15
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:15:02
% EndTime: 2019-02-22 10:15:02
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->15), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
t96 = cos(qJ(1));
t95 = sin(qJ(1));
t83 = sin(qJ(5));
t84 = sin(qJ(4));
t94 = t84 * t83;
t85 = cos(qJ(5));
t93 = t84 * t85;
t86 = cos(qJ(4));
t92 = t86 * t83;
t91 = t86 * t85;
t90 = cos(pkin(9));
t89 = sin(pkin(9));
t78 = -t95 * t89 - t96 * t90;
t79 = t96 * t89 - t95 * t90;
t88 = t78 * t83 + t79 * t91;
t87 = -t78 * t85 + t79 * t92;
t77 = -t78 * t91 + t79 * t83;
t76 = t78 * t92 + t79 * t85;
t1 = [t88, 0, 0, t78 * t93, t76, 0; t77, 0, 0, t79 * t93, t87, 0; 0, 0, 0, -t91, t94, 0; -t87, 0, 0, -t78 * t94, -t77, 0; t76, 0, 0, -t79 * t94, t88, 0; 0, 0, 0, t92, t93, 0; t79 * t84, 0, 0, -t78 * t86, 0, 0; -t78 * t84, 0, 0, -t79 * t86, 0, 0; 0, 0, 0, -t84, 0, 0;];
JR_rot  = t1;
